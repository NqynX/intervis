# intervis — Shiny GUI ---------------------------------------------------------
# The web tool: upload one genome for a detailed single-array view, or two to
# compare. It runs the intervis engine (IntegronFinder -> annotate -> cluster ->
# interactive viewer) and embeds the result. A thin front-end over the CLI, so
# the web tool and the command line stay in step.
#
# Launch (from the project root, inside the integronfinder env so integron_finder,
# diamond and python are on PATH):
#     micromamba activate integronfinder
#     R -e 'shiny::runApp("intervis/app", host="0.0.0.0", port=8787)'
# then open the forwarded port in your browser.
#
# Edit the three config lines below if your layout differs.

library(shiny)

# Shiny defaults to a 5 MB upload cap; genomes and GenBank files are bigger.
options(shiny.maxRequestSize = 500 * 1024^2)   # 500 MB

# runApp sets the wd to this app folder; the project root is one level up.
# Override any of these with env vars if your layout differs.
.app_dir     <- normalizePath(getwd())
PROJECT      <- Sys.getenv("INTERVIS_PROJECT",
                           unset = normalizePath(file.path(.app_dir, ".."), mustWork = FALSE))
INTERVIS_PKG <- Sys.getenv("INTERVIS_PKG",
                           unset = PROJECT)
SWISSPROT    <- Sys.getenv("INTERVIS_SWISSPROT",
                           unset = normalizePath(file.path(PROJECT, "db", "swissprot"), mustWork = FALSE))
CPU          <- as.integer(Sys.getenv("INTERVIS_CPU", unset = "8"))

# Bundled example datasets. AUTO-DISCOVERED from the examples/ folder — any
# <id>.fna (with optional <id>.gb/.gff and <id>.integrons) becomes an example, its
# menu label taken from the FASTA header. No manifest to be clobbered on update, so
# examples added by fetch_vibrio_examples.sh / add_example.sh persist across upgrades.
.ex_dir <- file.path(.app_dir, "examples")
.first <- function(cands) { hit <- cands[file.exists(file.path(.ex_dir, cands))]; if (length(hit)) hit[[1]] else "" }
.deflabel <- function(fasta, id) {
  if (!nzchar(fasta)) return(id)
  h <- tryCatch(readLines(file.path(.ex_dir, fasta), n = 1, warn = FALSE), error = function(e) "")
  if (!length(h) || !startsWith(h, ">")) return(id)
  acc  <- sub("^>(\\S*).*$", "\\1", h)
  desc <- sub("^>\\S*\\s*", "", h)
  desc <- trimws(sub(",?\\s*complete (sequence|genome|chromosome).*$", "", desc))
  # append the accession only when it looks like a real sequence accession (has a digit);
  # a synthetic demo id like "DemoVibrioA_chr" has none, so its label stays clean
  show_acc <- nzchar(acc) && grepl("[0-9]", acc)
  if (nzchar(desc)) paste0(desc, if (show_acc) paste0(" (", acc, ")") else "") else id
}
discover_examples <- function() {
  ff  <- list.files(.ex_dir, pattern = "\\.(fna|fasta|fa|integrons)$")
  ids <- unique(sub("\\.(fna|fasta|fa|integrons)$", "", ff)); ids <- sort(ids[nzchar(ids)])
  ex <- list()
  for (id in ids) {
    fasta <- .first(paste0(id, c(".fna", ".fasta", ".fa")))
    ex[[id]] <- list(id = id, label = .deflabel(fasta, id), fasta = fasta,
                     genbank   = .first(paste0(id, c(".gbff", ".gbk", ".gb", ".gff", ".gff3"))),
                     integrons = .first(paste0(id, ".integrons")))
  }
  ex
}
EX <- tryCatch(discover_examples(), error = function(e) list())
ex_path <- function(id, col) normalizePath(file.path(.ex_dir, EX[[id]][[col]]), mustWork = FALSE)
EX_CHOICES <- c("— none —" = "",
                if (length(EX)) setNames(vapply(EX, `[[`, "", "id"),
                                         vapply(EX, `[[`, "", "label")))

hlp <- function(txt)
  tags$p(txt, style = "color:#6b7280;font-size:11.5px;margin:2px 0 12px;line-height:1.45")

ui <- fluidPage(
  tags$head(tags$style(HTML(
    "body{font-family:system-ui,sans-serif} .well{background:#fbf9f6}
     #viewer iframe{width:100%;height:840px;border:1px solid #e6e1da;border-radius:12px}
     .run{background:#3b6ea5;color:#fff;font-weight:600;border:0}"))),
  titlePanel(title = div(
    span("intervis", style = "font-weight:700"),
    span("  ·  interactive integron array visualiser",
         style = "font-weight:400;color:#5b6b7c;font-size:16px")),
    windowTitle = "intervis"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$p(tags$b("Examples"), style = "margin-bottom:2px"),
      selectInput("ex1", "Example 1", choices = EX_CHOICES, selected = ""),
      selectInput("ex2", "Example 2 (optional, to compare)", choices = EX_CHOICES, selected = ""),
      hlp("Pick an example, or upload your own below."),
      tags$hr(),
      tags$p(tags$b("Genome (FASTA)"), style = "margin-bottom:2px"),
      fileInput("g1", "Genome 1", accept = c(".fna", ".fasta", ".fa")),
      fileInput("g2", "Genome 2 (optional)", accept = c(".fna", ".fasta", ".fa")),
      tags$hr(),
      tags$p(tags$b("Annotation (optional)"), style = "margin-bottom:2px"),
      fileInput("a1", "Genes 1 — GFF/GenBank", accept = c(".gff", ".gff3", ".gbff", ".gbk", ".gb")),
      fileInput("a2", "Genes 2 — GFF/GenBank", accept = c(".gff", ".gff3", ".gbff", ".gbk", ".gb")),
      tags$hr(),
      tags$p(tags$b("or .integrons"), style = "margin-bottom:2px"),
      fileInput("i1", "Array 1", accept = c(".integrons", ".tsv", ".txt")),
      fileInput("i2", "Array 2 (optional)", accept = c(".integrons", ".tsv", ".txt")),
      tags$hr(),
      checkboxInput("annotate", "Name cassettes (Swiss-Prot)", value = TRUE),
      checkboxInput("localmax", "Thorough attC search (--local-max)", value = TRUE),
      actionButton("run", "Run intervis", class = "run", width = "100%")
    ),
    mainPanel(width = 9, div(id = "viewer", uiOutput("viewer")))
  )
)

server <- function(input, output, session) {
  result <- eventReactive(input$run, {
    work <- tempfile("intervis_run_"); dir.create(work)
    out <- file.path(work, "result.html")
    old <- Sys.getenv("PYTHONPATH"); Sys.setenv(PYTHONPATH = INTERVIS_PKG)
    on.exit(Sys.setenv(PYTHONPATH = old), add = TRUE)
    cp <- function(fi) { p <- file.path(work, fi$name); file.copy(fi$datapath, p); p }

    ex1 <- if (!is.null(input$ex1)) input$ex1 else ""
    ex2 <- if (!is.null(input$ex2)) input$ex2 else ""

    if (!is.null(input$i1)) {                       # .integrons mode — instant, skips IntegronFinder
      i1 <- cp(input$i1)
      args <- if (!is.null(input$i2))
                c("-m", "intervis", "compare", i1, cp(input$i2), "-o", out)
              else
                c("-m", "intervis", "view", i1, "-o", out)
      detail <- "Rendering (instant)"
    } else if (is.null(input$g1) && nzchar(ex1)) {  # bundled example
      gbk     <- function(id) if (nzchar(EX[[id]]$genbank)) ex_path(id, "genbank") else "NONE"
      has_int <- function(id) nzchar(EX[[id]]$integrons) && file.exists(ex_path(id, "integrons"))
      cmp     <- nzchar(ex2)
      instant <- has_int(ex1) && (!cmp || has_int(ex2))
      if (instant && cmp) {                         # precomputed .integrons -> instant compare
        args <- c("-m", "intervis", "compare",
                  ex_path(ex1, "integrons"), ex_path(ex2, "integrons"),
                  "--genomes", ex_path(ex1, "fasta"), ex_path(ex2, "fasta"),
                  "--annotations", gbk(ex1), gbk(ex2), "-o", out, "--cpu", CPU)
        detail <- "Rendering example (instant)"
      } else if (instant) {                         # precomputed .integrons -> instant single view
        args <- c("-m", "intervis", "view", ex_path(ex1, "integrons"),
                  "--genome", ex_path(ex1, "fasta"), "-o", out)
        if (nzchar(EX[[ex1]]$genbank)) args <- c(args, "--annotation", ex_path(ex1, "genbank"))
        detail <- "Rendering example (instant)"
      } else {                                      # FASTA-only example -> run IntegronFinder
        genomes <- ex_path(ex1, "fasta"); if (cmp) genomes <- c(genomes, ex_path(ex2, "fasta"))
        args <- c("-m", "intervis", "run", genomes, "-o", out, "--outdir", work, "--cpu", CPU)
        an <- gbk(ex1); if (cmp) an <- c(an, gbk(ex2))
        if (any(an != "NONE")) args <- c(args, "--annotation", an)
        detail <- "IntegronFinder on the example (a few minutes) — watch the R terminal"
      }
    } else {                                        # genome mode — runs IntegronFinder
      req(input$g1)
      genomes <- cp(input$g1)
      if (!is.null(input$g2)) genomes <- c(genomes, cp(input$g2))
      args <- c("-m", "intervis", "run", genomes, "-o", out, "--outdir", work, "--cpu", CPU)
      if (input$annotate) args <- c(args, "--annotate", "--swissprot", SWISSPROT)
      if (!input$localmax) args <- c(args, "--no-local-max")
      # optional per-genome gene annotation (GFF/GenBank); "NONE" = none for that genome
      if (!is.null(input$a1) || !is.null(input$a2)) {
        an <- c(if (!is.null(input$a1)) cp(input$a1) else "NONE")
        if (!is.null(input$g2)) an <- c(an, if (!is.null(input$a2)) cp(input$a2) else "NONE")
        args <- c(args, "--annotation", an)
      }
      detail <- "IntegronFinder — watch the R terminal (this takes minutes)"
    }

    status <- 1L
    withProgress(message = "Running intervis", value = 0.4, {
      setProgress(0.5, detail = detail)
      message("\n[intervis] launching: python ", paste(args, collapse = " "))
      status <- system2("python", args, stdout = "", stderr = "")   # stream live to the R terminal
      setProgress(0.95, detail = "Rendering")
    })
    list(ok = file.exists(out),
         html = if (file.exists(out)) paste(readLines(out, warn = FALSE), collapse = "\n") else "",
         log = sprintf("intervis exited with status %s and produced no figure. The full log is in the R terminal where the app is running — the most common cause is a fragmented draft genome with --local-max on (untick 'Thorough attC search'), or IntegronFinder finding no integron in this genome.", status))
  })

  output$viewer <- renderUI({
    if (input$run == 0)
      return(div(style = "padding:40px 24px;color:#5b6b7c;max-width:560px",
        tags$h3("Pick an example, or upload, then press Run", style = "color:#1a2230;font-weight:600"),
        tags$p("New here? Choose a bundled Example and press Run to see a result instantly — "),
        tags$p("pick a second example to see the two-genome comparison. "),
        tags$p("Or give intervis your own genome (FASTA) and it runs IntegronFinder, names the "),
        tags$p("cassettes and draws the array; an existing .integrons file renders instantly.")))
    r <- result()
    if (!r$ok)
      return(div(style = "padding:24px;color:#b23a48",
        tags$h4("The run did not produce a figure."),
        tags$pre(style = "white-space:pre-wrap;font-size:12px;color:#5b6b7c", r$log)))
    tags$iframe(srcdoc = r$html)
  })
}

shinyApp(ui, server)

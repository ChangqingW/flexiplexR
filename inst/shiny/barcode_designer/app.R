library(shiny)
library(sortable)
library(jsonlite)

TYPE_COLORS <- c(
  FIXED         = "#9E9E9E",
  MATCHED       = "#42A5F5",
  RANDOM        = "#66BB6A",
  MATCHED_SPLIT = "#FFA726"
)

TYPE_DEFAULTS <- list(
  FIXED         = list(name = "primer",  pattern = "CTACACGACGCTCTTCCGATCT", bc_list_name = "", group = "", buffer_size = 2L, max_edit_distance = 2L),
  MATCHED       = list(name = "CB",      pattern = "NNNNNNNNNNNNNNNN",       bc_list_name = "CB", group = "", buffer_size = 5L, max_edit_distance = 2L),
  RANDOM        = list(name = "UB",      pattern = "NNNNNNNNNNNN",           bc_list_name = "", group = "", buffer_size = 2L, max_edit_distance = 2L),
  MATCHED_SPLIT = list(name = "BC1",     pattern = "NNNNNNNN",               bc_list_name = "", group = "CB", buffer_size = 2L, max_edit_distance = 2L)
)

# Card div — selection outline is applied by JS, not R
make_segment_card <- function(seg) {
  tags$div(
    class      = "segment-card",
    `data-uid` = as.character(seg$uid),
    style = paste0(
      "background:", TYPE_COLORS[[seg$type]], ";color:#fff;",
      "padding:8px 12px;margin:4px 0;border-radius:6px;",
      "cursor:grab;font-family:monospace;font-size:13px;",
      "display:flex;align-items:center;gap:10px;"
    ),
    tags$span(style = "font-weight:bold;min-width:90px;",              seg$name),
    tags$span(style = "opacity:0.85;min-width:70px;font-size:11px;",   paste0("[", seg$type, "]")),
    tags$span(style = "opacity:0.9;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;",
              seg$pattern)
  )
}

uid_to_idx <- function(segs, uid) {
  idx <- which(vapply(segs, `[[`, 0L, "uid") == uid)
  if (length(idx) == 0L) 1L else idx[[1L]]
}

segments_to_json <- function(segments, barcode_groups, params) {
  clean <- lapply(segments, function(s) { s$uid <- NULL; s })
  bp <- list(
    max_flank_editdistance  = params$max_flank_editdistance,
    strand                  = params$strand,
    TSO_seq                 = params$TSO_seq,
    TSO_prime               = params$TSO_prime,
    cutadapt_minimum_length = params$cutadapt_minimum_length,
    full_length_only        = params$full_length_only,
    segments                = clean,
    barcode_groups          = if (length(barcode_groups) == 0) list() else barcode_groups
  )
  jsonlite::toJSON(list(barcode_parameters = bp), pretty = TRUE, auto_unbox = TRUE)
}

subsample_fastq <- function(in_path, out_path, n_reads) {
  con <- if (grepl("\\.gz$", in_path, ignore.case = TRUE)) gzfile(in_path, "r") else file(in_path, "r")
  on.exit(close(con))
  lines <- readLines(con, n = 4L * n_reads)
  writeLines(lines, out_path)
}

segments_to_r <- function(segments, barcode_groups) {
  lines <- "segments <- list(\n"
  for (i in seq_along(segments)) {
    s     <- segments[[i]]
    comma <- if (i < length(segments)) "," else ""
    lines <- paste0(lines, switch(s$type,
      FIXED = sprintf(
        '  barcode_segment("FIXED", "%s", "%s")%s\n',
        s$pattern, s$name, comma),
      MATCHED = sprintf(
        '  barcode_segment("MATCHED", "%s", "%s", bc_list = "%s", buffer_size = %d, max_edit_distance = %d)%s\n',
        s$pattern, s$name, s$bc_list_name, s$buffer_size, s$max_edit_distance, comma),
      RANDOM = sprintf(
        '  barcode_segment("RANDOM", "%s", "%s")%s\n',
        s$pattern, s$name, comma),
      MATCHED_SPLIT = sprintf(
        '  barcode_segment("MATCHED_SPLIT", "%s", "%s", group = "%s", buffer_size = %d, max_edit_distance = %d)%s\n',
        s$pattern, s$name, s$group, s$buffer_size, s$max_edit_distance, comma)
    ))
  }
  lines <- paste0(lines, ")\n\n")
  if (length(barcode_groups) > 0) {
    lines <- paste0(lines, "barcode_groups <- list(\n")
    for (i in seq_along(barcode_groups)) {
      g     <- barcode_groups[[i]]
      comma <- if (i < length(barcode_groups)) "," else ""
      lines <- paste0(lines, sprintf(
        '  barcode_group(name = "%s", bc_list_name = "%s", max_edit_distance = %d)%s\n',
        g$name, g$bc_list_name, g$max_edit_distance, comma))
    }
    lines <- paste0(lines, ")\n")
  } else {
    lines <- paste0(lines, "barcode_groups <- list()\n")
  }
  lines
}

# ── UI ───────────────────────────────────────────────────────────────────────

.preview_fastq   <- shiny::getShinyOption("rflexiplex_fastq",          NULL)
.preview_bcs     <- shiny::getShinyOption("rflexiplex_barcodes_files", NULL)
.preview_enabled <- !is.null(.preview_fastq) && !is.null(.preview_bcs)

preview_tab <- tabPanel("Preview",
  br(),
  if (!.preview_enabled) {
    helpText(
      "Pass fastq and barcodes_files arguments to run_barcode_designer() ",
      "to enable a demultiplex preview here."
    )
  } else tagList(
    sliderInput("preview_n", "Reads to sample",
                min = 100, max = 10000, value = 1000, step = 100),
    actionButton("run_preview", "Run preview", class = "btn-primary"),
    br(), br(),
    uiOutput("preview_status"),
    tableOutput("preview_counts"),
    h5("Example demultiplexed read headers"),
    verbatimTextOutput("preview_headers"),
    h5("Example stats rows"),
    tableOutput("preview_stats")
  )
)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: sans-serif; }
      .segment-card { user-select: none; }
      .title-row { display:flex; align-items:center; justify-content:space-between; margin-bottom:16px; }
      .nav-tabs   { margin-bottom: 10px; }
    ")),
    tags$script(HTML("
      // Click a card → outline it and tell Shiny which uid is selected
      $(document).on('click', '.segment-card', function() {
        $('.segment-card').css('outline', '');
        $(this).css('outline', '3px solid white');
        var uid = parseInt($(this).attr('data-uid'));
        if (!isNaN(uid)) Shiny.setInputValue('clicked_uid', uid, {priority: 'event'});
      });

      // After renderUI re-renders the card list, re-apply the outline to the selected card
      Shiny.addCustomMessageHandler('highlightCard', function(uid) {
        $('.segment-card').css('outline', '');
        $('.segment-card[data-uid=\"' + uid + '\"]').css('outline', '3px solid white');
      });
    "))
  ),

  tags$div(class = "title-row",
    tags$h2(style = "margin:0;", "Barcode Segment Designer"),
    actionButton("done", "Done ✓", class = "btn-success btn-lg")
  ),

  fluidRow(
    # ── Left: card list + groups ─────────────────────────────────────────
    column(5,
      wellPanel(
        h4("Read Structure"),
        tags$p(style = "color:#555;font-size:12px;margin-bottom:6px;",
          "Drag to reorder. Click a card to select for editing."),
        # Cards rendered here; sortable_js is injected inside renderUI
        uiOutput("segment_cards"),
        br(),
        fluidRow(
          column(7, selectInput("new_type", NULL,
            choices = c("FIXED", "MATCHED", "RANDOM", "MATCHED_SPLIT"), width = "100%")),
          column(5, actionButton("add_segment", "Add", class = "btn-primary", width = "100%"))
        ),
        hr(),
        h4("Barcode Groups"),
        tags$p(style = "color:#555;font-size:12px;", "Required for MATCHED_SPLIT segments."),
        uiOutput("group_list_ui"),
        actionButton("add_group", "Add Group", class = "btn-default")
      )
    ),

    # ── Middle: editor ───────────────────────────────────────────────────
    column(3,
      wellPanel(
        h4("Edit Segment"),
        uiOutput("editor_ui"),
        hr(),
        h4("Global Parameters"),
        numericInput("max_flank_editdistance", "Max flank edit distance", value = 8, min = 0),
        selectInput("strand", "Strand", choices = c("+", "-"), selected = "-"),
        textInput("TSO_seq", "TSO sequence", value = "AAGCAGTGGTATCAACGCAGAGTACATGGG"),
        numericInput("TSO_prime", "TSO prime (3 or 5)", value = 5, min = 3, max = 5, step = 2),
        numericInput("cutadapt_minimum_length", "Min length after TSO trim", value = 10, min = 0),
        checkboxInput("full_length_only", "Full length only", value = FALSE)
      )
    ),

    # ── Right: preview ───────────────────────────────────────────────────
    column(4,
      wellPanel(
        tabsetPanel(
          tabPanel("JSON",
            br(),
            downloadButton("download_json", "Download JSON"),
            br(), br(),
            verbatimTextOutput("json_preview")
          ),
          tabPanel("R code",
            br(),
            verbatimTextOutput("r_preview")
          ),
          preview_tab
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Reactive state ────────────────────────────────────────────────────────
  segments <- reactiveVal(list(
    list(uid = 1L, type = "FIXED",   pattern = "CTACACGACGCTCTTCCGATCT", name = "primer",
         bc_list_name = "", group = "", buffer_size = 2L, max_edit_distance = 2L),
    list(uid = 2L, type = "MATCHED", pattern = "NNNNNNNNNNNNNNNN",       name = "CB",
         bc_list_name = "CB", group = "", buffer_size = 5L, max_edit_distance = 2L),
    list(uid = 3L, type = "RANDOM",  pattern = "NNNNNNNNNNNN",           name = "UB",
         bc_list_name = "", group = "", buffer_size = 2L, max_edit_distance = 2L),
    list(uid = 4L, type = "FIXED",   pattern = "TTTTTTTTT",              name = "polyT",
         bc_list_name = "", group = "", buffer_size = 2L, max_edit_distance = 2L)
  ))

  barcode_groups <- reactiveVal(list())
  next_uid       <- reactiveVal(5L)
  selected_uid   <- reactiveVal(1L)

  # ── Segment card list — depends ONLY on segments(), never selected_uid() ──
  #    Selection outline is managed by JavaScript to avoid re-render on click.
  output$segment_cards <- renderUI({
    segs <- segments()
    card_divs <- lapply(segs, make_segment_card)
    tagList(
      tags$div(id = "segment-card-list", card_divs),
      # Re-attach sortable every time this renderUI runs (after segments() changes)
      sortable_js(
        "segment-card-list",
        options = sortable_options(
          onEnd = htmlwidgets::JS("function(evt) {
            var items = document.getElementById('segment-card-list').children;
            var uids = Array.from(items).map(function(el) {
              return parseInt(el.getAttribute('data-uid'));
            });
            Shiny.setInputValue('segment_drag_order', uids, {priority: 'event'});
          }")
        )
      )
    )
  })

  # After a re-render, re-apply the selection outline via JS message
  observeEvent(segments(), {
    session$sendCustomMessage("highlightCard", selected_uid())
  }, ignoreInit = TRUE)

  # ── Drag-and-drop reorder ─────────────────────────────────────────────────
  observeEvent(input$segment_drag_order, {
    # JS sends integer UIDs in new DOM order via parseInt — no regex needed
    uids <- as.integer(input$segment_drag_order)
    segs <- segments()
    current_uids <- vapply(segs, `[[`, 0L, "uid")
    new_order    <- match(uids, current_uids)
    if (!anyNA(new_order) && length(new_order) == length(segs)) {
      segments(segs[new_order])
    }
  })

  # ── Click to select ───────────────────────────────────────────────────────
  observeEvent(input$clicked_uid, {
    selected_uid(input$clicked_uid)
  })

  # ── Add segment ───────────────────────────────────────────────────────────
  observeEvent(input$add_segment, {
    uid      <- next_uid()
    next_uid(uid + 1L)
    type     <- input$new_type
    new_seg  <- c(list(uid = uid, type = type), TYPE_DEFAULTS[[type]])
    segs     <- segments()
    segs[[length(segs) + 1L]] <- new_seg
    segments(segs)
    selected_uid(uid)
  })

  # ── Editor panel ──────────────────────────────────────────────────────────
  output$editor_ui <- renderUI({
    segs <- segments()
    if (length(segs) == 0L) return(p(style = "color:#888;", "No segments."))
    uid <- selected_uid()
    idx <- uid_to_idx(segs, uid)
    s   <- segs[[idx]]

    tagList(
      tags$div(
        style = paste0("padding:4px 8px;border-radius:4px;margin-bottom:10px;",
                       "background:", TYPE_COLORS[[s$type]], ";color:white;"),
        strong(paste0("[", idx, "/", length(segs), "]  ")),
        strong(s$type), " — ", s$name
      ),
      selectInput("edit_type", "Type",
        choices = c("FIXED", "MATCHED", "RANDOM", "MATCHED_SPLIT"), selected = s$type),
      textInput("edit_name",    "Name",    value = s$name),
      textInput("edit_pattern", "Pattern", value = s$pattern),
      conditionalPanel("input.edit_type == 'MATCHED'",
        textInput("edit_bc_list_name", "bc_list key (in barcodes_files)", value = s$bc_list_name)
      ),
      conditionalPanel("input.edit_type == 'MATCHED_SPLIT'",
        textInput("edit_group", "Group name", value = s$group)
      ),
      conditionalPanel("input.edit_type == 'MATCHED' || input.edit_type == 'MATCHED_SPLIT'",
        numericInput("edit_buffer_size",       "Buffer size",       value = s$buffer_size,       min = 0),
        numericInput("edit_max_edit_distance", "Max edit distance", value = s$max_edit_distance, min = 0)
      ),
      fluidRow(
        column(6, actionButton("apply_edit",     "Apply",  class = "btn-success", width = "100%")),
        column(6, actionButton("delete_segment", "Delete", class = "btn-danger",  width = "100%"))
      )
    )
  })

  observeEvent(input$apply_edit, {
    segs <- segments()
    uid  <- selected_uid()
    idx  <- uid_to_idx(segs, uid)
    type <- input$edit_type
    segs[[idx]] <- list(
      uid               = uid,
      type              = type,
      name              = input$edit_name,
      pattern           = input$edit_pattern,
      bc_list_name      = if (type == "MATCHED")       input$edit_bc_list_name
                          else if (type == "MATCHED_SPLIT") input$edit_group
                          else "",
      group             = if (type == "MATCHED_SPLIT") input$edit_group else "",
      buffer_size       = if (type %in% c("MATCHED", "MATCHED_SPLIT")) as.integer(input$edit_buffer_size)       else 2L,
      max_edit_distance = if (type %in% c("MATCHED", "MATCHED_SPLIT")) as.integer(input$edit_max_edit_distance) else 2L
    )
    segments(segs)
  })

  observeEvent(input$delete_segment, {
    segs <- segments()
    uid  <- selected_uid()
    idx  <- uid_to_idx(segs, uid)
    segs[[idx]] <- NULL
    segments(segs)
    if (length(segs) > 0L) {
      selected_uid(segs[[min(idx, length(segs))]]$uid)
    }
  })

  # ── Barcode groups ────────────────────────────────────────────────────────
  output$group_list_ui <- renderUI({
    grps <- barcode_groups()
    if (length(grps) == 0L) return(p(style = "color:#888;", "No barcode groups."))
    lapply(seq_along(grps), function(i) {
      g <- grps[[i]]
      tags$div(style = "display:flex;align-items:center;gap:8px;margin:4px 0;",
        tags$span(
          style = "background:#FFA726;color:white;padding:3px 8px;border-radius:4px;font-size:12px;",
          paste0(g$name, " → ", g$bc_list_name, " (ed≤", g$max_edit_distance, ")")
        ),
        actionButton(paste0("del_group_", i), "×", class = "btn-xs btn-danger")
      )
    })
  })

  observeEvent(input$add_group, {
    showModal(modalDialog(
      title = "Add Barcode Group",
      textInput("new_group_name",    "Group name",         value = "CB"),
      textInput("new_group_bc_list", "bc_list_name",       value = "CB"),
      numericInput("new_group_max_ed", "Max edit distance", value = 2, min = 0),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_add_group", "Add", class = "btn-primary")
      )
    ))
  })

  observeEvent(input$confirm_add_group, {
    grps <- barcode_groups()
    grps[[length(grps) + 1L]] <- list(
      name              = input$new_group_name,
      bc_list_name      = input$new_group_bc_list,
      max_edit_distance = as.integer(input$new_group_max_ed)
    )
    barcode_groups(grps)
    removeModal()
  })

  observe({
    grps <- barcode_groups()
    lapply(seq_along(grps), function(i) {
      observeEvent(input[[paste0("del_group_", i)]], {
        g2      <- barcode_groups()
        g2[[i]] <- NULL
        barcode_groups(g2)
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # ── Global params ─────────────────────────────────────────────────────────
  global_params <- reactive({
    list(
      max_flank_editdistance  = input$max_flank_editdistance,
      strand                  = input$strand,
      TSO_seq                 = input$TSO_seq,
      TSO_prime               = input$TSO_prime,
      cutadapt_minimum_length = input$cutadapt_minimum_length,
      full_length_only        = input$full_length_only
    )
  })

  # ── Previews ──────────────────────────────────────────────────────────────
  json_text <- reactive({
    segments_to_json(segments(), barcode_groups(), global_params())
  })

  output$json_preview  <- renderText({ json_text() })

  output$download_json <- downloadHandler(
    filename = "barcode_config.json",
    content  = function(f) writeLines(json_text(), f)
  )

  output$r_preview <- renderText({
    segments_to_r(segments(), barcode_groups())
  })

  # ── Preview ───────────────────────────────────────────────────────────────
  if (.preview_enabled) {
    observeEvent(input$run_preview, {
      output$preview_status <- renderUI(
        div(style = "color:#555;", em("Running flexiplex on subsample..."))
      )
      tmp_fq    <- tempfile(fileext = ".fastq")
      tmp_out   <- tempfile(fileext = ".fastq")
      tmp_stats <- tempfile(fileext = ".tsv")
      session$onSessionEnded(function() {
        try(unlink(c(tmp_fq, tmp_out, tmp_stats)), silent = TRUE)
      })

      res <- tryCatch({
        subsample_fastq(.preview_fastq, tmp_fq, input$preview_n)
        seg_objs <- lapply(segments(), function(s) {
          Rflexiplex::barcode_segment(
            type    = s$type,
            pattern = s$pattern,
            name    = s$name,
            bc_list = if (isTRUE(nzchar(s$bc_list_name))) s$bc_list_name else NA_character_,
            group   = if (isTRUE(nzchar(s$group)))        s$group        else NA_character_,
            buffer_size       = s$buffer_size,
            max_edit_distance = s$max_edit_distance
          )
        })
        grp_objs <- lapply(barcode_groups(), function(g) {
          Rflexiplex::barcode_group(
            name              = g$name,
            bc_list_name      = g$bc_list_name,
            max_edit_distance = g$max_edit_distance
          )
        })
        Rflexiplex::find_barcode(
          fastq                   = tmp_fq,
          segments                = seg_objs,
          barcode_groups          = grp_objs,
          barcodes_files          = .preview_bcs,
          max_flank_editdistance  = input$max_flank_editdistance,
          reads_out               = tmp_out,
          stats_out               = tmp_stats,
          threads                 = 1,
          TSO_seq                 = if (is.null(input$TSO_seq)) "" else input$TSO_seq,
          TSO_prime               = input$TSO_prime,
          strand                  = input$strand,
          cutadapt_minimum_length = input$cutadapt_minimum_length,
          full_length_only        = input$full_length_only
        )
      }, error = function(e) e)

      if (inherits(res, "error")) {
        output$preview_status <- renderUI(
          div(style = "color:#c00;", strong("Error: "), res$message)
        )
        output$preview_counts  <- renderTable(NULL)
        output$preview_headers <- renderText("")
        output$preview_stats   <- renderTable(NULL)
        return()
      }

      rc    <- res$read_counts
      total <- as.integer(rc[["total reads"]])
      demux <- as.integer(rc[["demultiplexed reads"]])
      rate  <- if (total > 0) round(100 * demux / total, 2) else 0
      output$preview_status <- renderUI(
        div(style = "color:#2e7d32;", strong("Done. "),
            sprintf("%d/%d reads demultiplexed (%.2f%%).", demux, total, rate))
      )
      output$preview_counts <- renderTable(
        data.frame(metric = c(names(rc), "demux rate (%)"),
                   value  = c(as.integer(rc), rate)),
        striped = TRUE, rownames = FALSE
      )

      hdr_lines <- tryCatch(readLines(tmp_out, n = 20), error = function(e) character(0))
      headers   <- if (length(hdr_lines) >= 1)
        hdr_lines[seq(1, length(hdr_lines), by = 4)] else character(0)
      output$preview_headers <- renderText(
        if (length(headers)) paste(headers, collapse = "\n") else "(no demultiplexed reads)"
      )

      stats_df <- tryCatch(utils::head(read.delim(tmp_stats), 10),
                           error = function(e) data.frame())
      output$preview_stats <- renderTable(stats_df, striped = TRUE, rownames = FALSE)
    })
  }

  # ── Done ──────────────────────────────────────────────────────────────────
  observeEvent(input$done, {
    clean_segs <- lapply(segments(), function(s) { s$uid <- NULL; s })
    shiny::stopApp(returnValue = list(
      segments       = clean_segs,
      barcode_groups = barcode_groups()
    ))
  })
}

shinyApp(ui, server)

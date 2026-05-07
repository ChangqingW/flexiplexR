library(shiny)
library(sortable)
library(jsonlite)

# Segment type colours
TYPE_COLORS <- c(
  FIXED         = "#9E9E9E",
  MATCHED       = "#42A5F5",
  RANDOM        = "#66BB6A",
  MATCHED_SPLIT = "#FFA726"
)

# Default field values by type
TYPE_DEFAULTS <- list(
  FIXED         = list(name = "primer",  pattern = "CTACACGACGCTCTTCCGATCT", bc_list_name = "", group = "", buffer_size = 2, max_edit_distance = 2),
  MATCHED       = list(name = "CB",      pattern = "NNNNNNNNNNNNNNNN",       bc_list_name = "CB", group = "", buffer_size = 5, max_edit_distance = 2),
  RANDOM        = list(name = "UB",      pattern = "NNNNNNNNNNNN",           bc_list_name = "", group = "", buffer_size = 2, max_edit_distance = 2),
  MATCHED_SPLIT = list(name = "BC1",     pattern = "NNNNNNNN",               bc_list_name = "", group = "CB", buffer_size = 2, max_edit_distance = 2)
)

make_segment_card <- function(seg, idx) {
  col <- TYPE_COLORS[[seg$type]]
  tags$div(
    class = "segment-card",
    style = paste0(
      "background:", col, ";color:#fff;padding:8px 12px;margin:4px 0;",
      "border-radius:6px;cursor:grab;font-family:monospace;font-size:13px;",
      "display:flex;align-items:center;gap:10px;"
    ),
    tags$span(style = "font-weight:bold;min-width:100px;", seg$name),
    tags$span(style = "opacity:0.85;min-width:60px;font-size:11px;", paste0("[", seg$type, "]")),
    tags$span(style = "opacity:0.9;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;", seg$pattern)
  )
}

segments_to_json <- function(segments, barcode_groups, params) {
  bp <- list(
    max_flank_editdistance = params$max_flank_editdistance,
    strand                 = params$strand,
    TSO_seq                = params$TSO_seq,
    TSO_prime              = params$TSO_prime,
    cutadapt_minimum_length = params$cutadapt_minimum_length,
    full_length_only       = params$full_length_only,
    segments               = segments,
    barcode_groups         = if (length(barcode_groups) == 0) list() else barcode_groups
  )
  jsonlite::toJSON(list(barcode_parameters = bp), pretty = TRUE, auto_unbox = TRUE)
}

segments_to_r <- function(segments, barcode_groups) {
  lines <- "segments <- list(\n"
  for (i in seq_along(segments)) {
    s <- segments[[i]]
    comma <- if (i < length(segments)) "," else ""
    if (s$type == "FIXED") {
      lines <- paste0(lines, sprintf(
        '  barcode_segment("FIXED", "%s", "%s")%s\n',
        s$pattern, s$name, comma
      ))
    } else if (s$type == "MATCHED") {
      lines <- paste0(lines, sprintf(
        '  barcode_segment("MATCHED", "%s", "%s", bc_list = "%s", buffer_size = %d, max_edit_distance = %d)%s\n',
        s$pattern, s$name, s$bc_list_name, s$buffer_size, s$max_edit_distance, comma
      ))
    } else if (s$type == "RANDOM") {
      lines <- paste0(lines, sprintf(
        '  barcode_segment("RANDOM", "%s", "%s")%s\n',
        s$pattern, s$name, comma
      ))
    } else if (s$type == "MATCHED_SPLIT") {
      lines <- paste0(lines, sprintf(
        '  barcode_segment("MATCHED_SPLIT", "%s", "%s", group = "%s", buffer_size = %d, max_edit_distance = %d)%s\n',
        s$pattern, s$name, s$group, s$buffer_size, s$max_edit_distance, comma
      ))
    }
  }
  lines <- paste0(lines, ")\n\n")

  if (length(barcode_groups) > 0) {
    lines <- paste0(lines, "barcode_groups <- list(\n")
    for (i in seq_along(barcode_groups)) {
      g <- barcode_groups[[i]]
      comma <- if (i < length(barcode_groups)) "," else ""
      lines <- paste0(lines, sprintf(
        '  barcode_group(name = "%s", bc_list_name = "%s", max_edit_distance = %d)%s\n',
        g$name, g$bc_list_name, g$max_edit_distance, comma
      ))
    }
    lines <- paste0(lines, ")\n")
  } else {
    lines <- paste0(lines, "barcode_groups <- list()\n")
  }
  lines
}

# ── UI ──────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: sans-serif; }
    .segment-card { user-select: none; }
    .panel-heading { font-weight: bold; }
    #segment_order .rank-list-item { padding: 0 !important; background: transparent !important; border: none !important; }
    .nav-tabs { margin-bottom: 10px; }
  "))),
  titlePanel("Barcode Segment Designer"),
  fluidRow(
    # ── Left: Segment list (drag-and-drop) ────────────────────────────────
    column(5,
      wellPanel(
        h4("Read Structure"),
        p(style = "color:#555;font-size:12px;", "Drag to reorder segments. Click a segment to edit."),
        uiOutput("segment_rank_list"),
        fluidRow(
          column(6, selectInput("new_type", NULL,
            choices = c("FIXED", "MATCHED", "RANDOM", "MATCHED_SPLIT"), width = "100%")),
          column(6, actionButton("add_segment", "Add Segment", class = "btn-primary", width = "100%"))
        ),
        hr(),
        h4("Barcode Groups"),
        p(style = "color:#555;font-size:12px;", "Required when using MATCHED_SPLIT segments."),
        uiOutput("group_list_ui"),
        actionButton("add_group", "Add Group", class = "btn-default")
      )
    ),

    # ── Middle: Segment editor ────────────────────────────────────────────
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

    # ── Right: Output preview ─────────────────────────────────────────────
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
          )
        )
      )
    )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Reactive state
  segments <- reactiveVal(list(
    list(type = "FIXED",   pattern = "CTACACGACGCTCTTCCGATCT", name = "primer",
         bc_list_name = "", group = "", buffer_size = 2, max_edit_distance = 2),
    list(type = "MATCHED", pattern = "NNNNNNNNNNNNNNNN",       name = "CB",
         bc_list_name = "CB", group = "", buffer_size = 5, max_edit_distance = 2),
    list(type = "RANDOM",  pattern = "NNNNNNNNNNNN",           name = "UB",
         bc_list_name = "", group = "", buffer_size = 2, max_edit_distance = 2),
    list(type = "FIXED",   pattern = "TTTTTTTTT",              name = "polyT",
         bc_list_name = "", group = "", buffer_size = 2, max_edit_distance = 2)
  ))

  barcode_groups <- reactiveVal(list())
  selected_idx   <- reactiveVal(NULL)

  # ── Segment rank list ────────────────────────────────────────────────────
  output$segment_rank_list <- renderUI({
    segs <- segments()
    items <- lapply(seq_along(segs), function(i) make_segment_card(segs[[i]], i))
    rank_list(
      text       = NULL,
      labels     = items,
      input_id   = "segment_order",
      css_id     = "segment_order_list",
      class      = "default-sortable"
    )
  })

  # Reorder when drag-and-drop fires
  observeEvent(input$segment_order, {
    order_labels <- input$segment_order
    if (is.null(order_labels)) return()
    segs <- segments()
    # rank_list returns labels as character; we identify positions by card content
    # Instead, we track order via the label index
    new_order <- as.integer(order_labels)
    if (!anyNA(new_order) && length(new_order) == length(segs)) {
      segments(segs[new_order])
    }
  }, ignoreInit = TRUE)

  # ── Add segment ──────────────────────────────────────────────────────────
  observeEvent(input$add_segment, {
    type <- input$new_type
    defaults <- TYPE_DEFAULTS[[type]]
    new_seg <- c(list(type = type), defaults)
    segs <- segments()
    segs[[length(segs) + 1]] <- new_seg
    segments(segs)
    selected_idx(length(segs))
  })

  # ── Editor: click on a segment to edit ───────────────────────────────────
  # Because rank_list doesn't fire click events directly, we use numeric input for selection
  output$editor_ui <- renderUI({
    segs <- segments()
    if (length(segs) == 0) return(p("No segments."))

    idx <- selected_idx()
    if (is.null(idx) || idx > length(segs)) idx <- 1
    s <- segs[[idx]]

    tagList(
      numericInput("edit_idx", "Select segment (index)", value = idx, min = 1, max = length(segs), step = 1),
      tags$div(style = paste0("padding:4px 8px;border-radius:4px;background:", TYPE_COLORS[[s$type]], ";color:white;margin-bottom:8px;"),
        strong(s$type), " — ", s$name),
      selectInput("edit_type", "Type", choices = c("FIXED", "MATCHED", "RANDOM", "MATCHED_SPLIT"), selected = s$type),
      textInput("edit_name", "Name", value = s$name),
      textInput("edit_pattern", "Pattern", value = s$pattern),
      conditionalPanel("input.edit_type == 'MATCHED'",
        textInput("edit_bc_list_name", "bc_list (key in barcodes_files)", value = s$bc_list_name)
      ),
      conditionalPanel("input.edit_type == 'MATCHED_SPLIT'",
        textInput("edit_group", "Group name", value = s$group)
      ),
      conditionalPanel("input.edit_type == 'MATCHED' || input.edit_type == 'MATCHED_SPLIT'",
        numericInput("edit_buffer_size", "Buffer size", value = s$buffer_size, min = 0),
        numericInput("edit_max_edit_distance", "Max edit distance", value = s$max_edit_distance, min = 0)
      ),
      fluidRow(
        column(6, actionButton("apply_edit", "Apply", class = "btn-success", width = "100%")),
        column(6, actionButton("delete_segment", "Delete", class = "btn-danger", width = "100%"))
      )
    )
  })

  observeEvent(input$edit_idx, {
    selected_idx(as.integer(input$edit_idx))
  })

  observeEvent(input$apply_edit, {
    idx <- as.integer(input$edit_idx)
    segs <- segments()
    if (idx < 1 || idx > length(segs)) return()
    type <- input$edit_type
    segs[[idx]] <- list(
      type             = type,
      name             = input$edit_name,
      pattern          = input$edit_pattern,
      bc_list_name     = if (type == "MATCHED") input$edit_bc_list_name else if (type == "MATCHED_SPLIT") input$edit_group else "",
      group            = if (type == "MATCHED_SPLIT") input$edit_group else "",
      buffer_size      = if (type %in% c("MATCHED", "MATCHED_SPLIT")) as.integer(input$edit_buffer_size) else 2L,
      max_edit_distance = if (type %in% c("MATCHED", "MATCHED_SPLIT")) as.integer(input$edit_max_edit_distance) else 2L
    )
    segments(segs)
  })

  observeEvent(input$delete_segment, {
    idx <- as.integer(input$edit_idx)
    segs <- segments()
    if (length(segs) == 0) return()
    segs[[idx]] <- NULL
    segments(segs)
    selected_idx(min(idx, length(segs)))
  })

  # ── Barcode groups ───────────────────────────────────────────────────────
  output$group_list_ui <- renderUI({
    grps <- barcode_groups()
    if (length(grps) == 0) return(p(style = "color:#888;", "No barcode groups."))
    lapply(seq_along(grps), function(i) {
      g <- grps[[i]]
      tags$div(style = "display:flex;align-items:center;gap:8px;margin:4px 0;",
        tags$span(style = "background:#FFA726;color:white;padding:3px 8px;border-radius:4px;font-size:12px;",
          paste0(g$name, " → ", g$bc_list_name)),
        actionButton(paste0("del_group_", i), "×", class = "btn-xs btn-danger")
      )
    })
  })

  observeEvent(input$add_group, {
    showModal(modalDialog(
      title = "Add Barcode Group",
      textInput("new_group_name", "Group name", value = "CB"),
      textInput("new_group_bc_list", "bc_list_name", value = "CB"),
      numericInput("new_group_max_ed", "Max edit distance", value = 2, min = 0),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_add_group", "Add", class = "btn-primary")
      )
    ))
  })

  observeEvent(input$confirm_add_group, {
    grps <- barcode_groups()
    grps[[length(grps) + 1]] <- list(
      name              = input$new_group_name,
      bc_list_name      = input$new_group_bc_list,
      max_edit_distance = as.integer(input$new_group_max_ed)
    )
    barcode_groups(grps)
    removeModal()
  })

  # Delete group buttons (dynamic)
  observe({
    grps <- barcode_groups()
    lapply(seq_along(grps), function(i) {
      btn_id <- paste0("del_group_", i)
      observeEvent(input[[btn_id]], {
        g2 <- barcode_groups()
        g2[[i]] <- NULL
        barcode_groups(g2)
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # ── Global params helper ─────────────────────────────────────────────────
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

  # ── JSON preview ─────────────────────────────────────────────────────────
  json_text <- reactive({
    segments_to_json(segments(), barcode_groups(), global_params())
  })

  output$json_preview <- renderText({ json_text() })

  output$download_json <- downloadHandler(
    filename = "barcode_config.json",
    content  = function(f) writeLines(json_text(), f)
  )

  # ── R code preview ───────────────────────────────────────────────────────
  output$r_preview <- renderText({
    segments_to_r(segments(), barcode_groups())
  })
}

shinyApp(ui, server)

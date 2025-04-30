suppressPackageStartupMessages({
  library(Seurat)
  library(purrr)
  library(forcats)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  library(enrichplot)
  library(CellChat)
  library(slingshot)
})


ptp <- function(
    pseudo.out,
    embedding,
    show.curves = TRUE,
    subset.expr = TRUE,
    path.arrow = arrow(),
    label.by = pseudo.out[["group.by"]],
    non.lineage.facet = NULL,
    line.palette = "rocket",
    point.palette = "turbo") {
  cell.dt <- pseudo.out[["input.embeddings"]][
    , .(
      cell,
      dim_1 = get(sprintf("%s_1", embedding)),
      dim_2 = get(sprintf("%s_2", embedding))
    )
  ][
    pseudo.out[["input.meta"]],
    on = .(cell)
  ]
  pseudo.dt <- pseudo.out[["pseudotime"]][
    cell.dt,
    on = .(cell)
  ][eval(substitute(subset.expr))][!is.na(pseudotime)]

  if (!is.null(non.lineage.facet) && uniqueN(pseudo.dt, "lineage") > 1) {
    stop("cannot have more than one lineage present when specifying a non-lineage faceting variable")
  }

  ggplot() +
    aes(x = dim_1, y = dim_2) +
    geom_point(
      aes(fill = pseudotime),
      pseudo.dt,
      stroke = 0,
      shape = "circle filled"
    ) +
    (if (show.curves) {
      geom_path(
        aes(colour = ord, group = NULL),
        pseudo.out[[
          sprintf("curves.%s", embedding)
        ]][
          lineage %in% pseudo.dt[, unique(lineage)]
        ],
        linewidth = 1,
        arrow = path.arrow
      )
    }) +
    (if (!is.null(label.by)) {
      geom_text(
        aes(label = cc),
        cell.dt[
          !is.na(get(label.by)),
          as.data.frame.list(colMeans(.SD[, .(dim_1, dim_2)])),
          by = .(cc = get(label.by))
        ]
      )
    }) +
    facet_wrap(
      vars(get(if (is.null(non.lineage.facet)) {
        "lineage"
      } else {
        non.lineage.facet
      })),
      axes = "all",
      axis.labels = "all"
    ) +
    scale_fill_viridis_c(option = point.palette) +
    scale_colour_viridis_c(option = line.palette) +
    labs(fill = "pseudotime", colour = "order")
}

mwb <- WhiteBackground(
  legend.background = element_rect(fill = "white"),
  legend.text = element_text(colour = "black"),
  legend.title = element_text(colour = "black"),
  legend.box.background = element_blank()
)

fig.1.c <- function(so) {
  sprintf("%s %s", so[[]][["genotype"]], so[[]][["condition"]]) |>
    fct_relevel("WT uninf.", "PINK1 KO uninf.", "WT inf.", "PINK1 KO inf.") |>
    AddMetaData(so, metadata = _, "split") |>
    DimPlot(split.by = "split", group.by = "timepoint")
}

fig.1.d <- function(so) {
  BC <- c("Ms4a1", "Cd19", "Pax5", "Bank1")
  DC <- c("Cd86", "Itgax", "Flt3", "Wdfy4")
  Mye <- c("Itgam", "Csf1r", "Adgre1", "Mrc1")
  NK <- c("Kit", "Klrb1b", "Klrk1")
  Gran <- c("Hdc", "Fpr2")
  gdT <- c("Gata3", "Tcrg-V1")
  CD4 <- c(
    "Cd3d", "Cd3e", "Cd4", "Tbx21",
    "Rorc", "Foxp3", "Ifng", "Il4", "Il10", "Il17a"
  )
  CD8 <- c("Cd8a", "Cd8b1", "Eomes", "Prf1", "Gzmb")
  naive <- c("Sell", "Ccr7", "Cd69", "Tcf7")
  resident <- c("Itgae")
  clusters_oi <- c(
    "B cell", "CD4 Th17", "Monocyte/M\u0278",
    "CD4 Th1", "CD8 Tc1", "CD8 Tn", "Granulocyte",
    "CD4 Th2", "gdT cell", "Dendritic cell",
    "NK/ILC", "dnT cell", "CD4 Treg",
    "CD4 Tn", "CD4 Trm", "CD8 Trm"
  )

  DotPlot(
    subset(so, celltype %in% clusters_oi),
    c(BC, DC, Mye, NK, Gran, gdT, CD4, CD8, naive, resident),
    group.by = "celltype", col.min = 0
  ) + scale_colour_viridis_c(option = "cividis") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    labs(y = NULL)
}

fig.1.e <- function(so) {
  DT <- as.data.table(Embeddings(so, "umap"), keep.rownames = "cell")[
    as.data.table(so[[c("timepoint", "celltype")]], keep.rownames = "cell"),
    on = .(cell)
  ]
  ggplot() +
    aes(x = umap_1, y = umap_2) +
    geom_point(aes(colour = timepoint), DT) +
    geom_text(aes(label = celltype), DT[, .(
      umap_1 = mean(umap_1),
      umap_2 = mean(umap_2)
    ), by = .(celltype)]) +
    mwb +
    theme(axis.ticks = element_blank(), axis.text = element_blank()) +
    labs(x = NULL, y = NULL)
}

fig.2.a <- function(so) {
  DT <- as.data.table(Embeddings(so, "umap"), keep.rownames = "cell")[
    as.data.table(so[[c(
      "timepoint", "genotype",
      "condition", "celltype"
    )]], keep.rownames = "cell"),
    on = .(cell)
  ][, `:=`(
    group = fifelse(condition == "uninf.", condition, timepoint)
  )][]
  ggplot() +
    aes(x = umap_1, y = umap_2) +
    geom_point(aes(colour = celltype), DT, show.legend = FALSE) +
    geom_text(aes(label = celltype), DT[, .(
      umap_1 = mean(umap_1),
      umap_2 = mean(umap_2)
    ), by = .(celltype)]) +
    mwb +
    theme(axis.ticks = element_blank(), axis.text = element_blank()) +
    labs(x = NULL, y = NULL) +
    facet_wrap(vars(group, genotype), ncol = 2, nrow = 3)
}

fig.2.b <- function(so) {
  clusters_oi <- c(
    "Monocyte/M\u0278", "Dendritic cell", "Granulocyte", "NK/ILC", "B cell",
    "gdT cell", "dnT cell", "CD4 Tn", "CD4 Trm", "CD4 Th1",
    "CD4 Th2", "CD4 Th17", "CD8 Tc1", "CD8 Tn", "CD8 Trm", "CD4 Treg"
  )
  list(
    "uninf." = subset(so, condition == "uninf."),
    "1-wpi" = subset(so, condition == "inf." & timepoint == "1-wpi"),
    "2-wpi" = subset(so, condition == "inf." & timepoint == "2-wpi")
  ) |>
    imap(function(so, group) {
      map(clusters_oi, function(cluster) {
        tryCatch(
          FindMarkers(
            subset(so, celltype == cluster),
            ident.1 = "PINK1 KO",
            ident.2 = "WT",
            group.by = "genotype",
            log2fc.threshold = 0.25
          ),
          error = function(error.cond) {
            print(sprintf("no DEGs for: %s %s", cluster, group))
            data.table(
              p_val = numeric(),
              avg_log2FC = numeric(),
              pct.1 = numeric(),
              pct.2 = numeric(),
              p_val_adj = numeric()
            )
          }
        ) |>
          as.data.table() |>
          _[
            p_val_adj < 0.05,
            .(sign = fifelse(avg_log2FC > 0, "pos", "neg"))
          ][
            ,
            .(n = .N),
            by = .(sign)
          ][
            data.table(n = 0, sign = c("pos", "neg")),
            .(
              n = fcoalesce(i.n + n, 0),
              timepoint = group,
              celltype = cluster,
              sign
            ),
            on = .(sign)
          ]
      }) |>
        reduce(rbind)
    }) |>
    reduce(rbind) |>
    dcast(timepoint + celltype ~ sign, value.var = "n") |>
    _[order(fct_relevel(timepoint, "uninf.", "1-wpi", "2-wpi"))] |>
    split(by = "timepoint") |>
    imap(function(data, timepoint) {
      p <- function(metric, colour) {
        ggplot(data) +
          aes(
            x = do.call(partial(fct_relevel, celltype), as.list(clusters_oi)),
            y = NA
          ) +
          geom_tile(aes(fill = get(metric))) +
          geom_text(aes(label = get(metric))) +
          labs(y = NULL) +
          scale_fill_gradient(
            low = "white",
            high = colour,
            breaks = c(100, 300),
            limits = c(0, 400),
            name = sprintf("%s-vly regulated", metric)
          )
      }
      list(p("pos", "red"), p("neg", "blue"))
    }) |>
    list_flatten() |>
    wrap_plots(... = _, guides = "collect", nrow = 6) +
    plot_layout(axes = "collect_x", axis_titles = "collect_x") &
    labs(x = NULL) &
    theme_void() &
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90, margin = margin(t = 2)),
      plot.background = element_rect(fill = "white"),
      axis.line.x = element_line(colour = "black")
    )
}

fig.3.a <- function(cd4.pt, cd8.pt, plot.embedding = "pca") {
  wrap_plots(
    ptp(cd4.pt, plot.embedding, subset.expr = lineage == "3") +
      guides(colour = "none") +
      theme_void() +
      theme(strip.text = element_blank()),
    ptp(cd4.pt, plot.embedding,
      subset.expr = lineage == "3",
      non.lineage.facet = "genotype",
      show.curves = FALSE, label.by = NULL
    ) + theme_void(),
    ptp(cd8.pt, plot.embedding) +
      guides(colour = "none") +
      theme_void() +
      theme(strip.text = element_blank()),
    ptp(cd8.pt, plot.embedding,
      non.lineage.facet = "genotype",
      show.curves = FALSE, label.by = NULL
    ) +
      theme_void(),
    nrow = 2,
    widths = c(1, 2)
  ) & guides(colour = "none")
}

mk_ratio_fig <- function(cd.stats, eff_type) {
  cd.stats[
    celltype.fine == "Tn"
  ][
    cd.stats[celltype.fine == eff_type],
    .(ratio = i.n / n, genotype = genotype),
    on = .(genotype)
  ] |>
    ggplot() +
    aes(x = fct_relevel(genotype, "WT"), y = ratio, fill = genotype) +
    geom_col() +
    labs(title = eff_type, y = "Teff/Tn ratio", x = NULL) +
    guides(fill = "none")
}

# only infected cells 2-wpi
fig.3.b <- function(cd4.inf.so, cd8.inf.so) {
  cd4.stats <- cd4.inf.so[[c("genotype", "celltype.fine")]] |>
    as.data.table() |>
    _[, .(n = .N), by = .(celltype.fine, genotype)]
  cd8.stats <- cd8.inf.so[[c("genotype", "celltype.fine")]] |>
    as.data.table() |>
    _[, .(n = .N), by = .(celltype.fine, genotype)]

  wrap_plots(
    ggplot(cd4.stats) +
      aes(x = fct_relevel(genotype, "WT"), y = n, fill = celltype.fine) +
      geom_col(position = position_fill()) +
      guides(fill = guide_legend(title = "CD4")) +
      theme_void() +
      theme(axis.text.x = element_text()),
    mk_ratio_fig(cd4.stats, "Th1"),
    mk_ratio_fig(cd4.stats, "Th17"),
    ggplot(cd8.stats) +
      aes(x = fct_relevel(genotype, "WT"), y = n, fill = celltype.fine) +
      geom_col(position = position_fill()) +
      guides(fill = guide_legend(title = "CD4")) +
      theme_void() +
      theme(axis.text.x = element_text()),
    mk_ratio_fig(cd8.stats, "Tc1"),
    mk_ratio_fig(cd8.stats, "Trm")
  )
}

fig.3.c <- function(th1.eGO) {
  wrap_plots(
    treeplot(th1.eGO, showCategory = 40),
    cnetplot(th1.eGO,
      showCategory = 40,
      max.overlaps = Inf, node_label = "gene"
    )
  )
}

fig.3.d <- function(tc1.eGO) {
  wrap_plots(
    treeplot(tc1.eGO, showCategory = 40),
    cnetplot(tc1.eGO,
      showCategory = 40,
      max.overlaps = Inf, node_label = "gene"
    )
  )
}

fig.4.a <- function(so.inf) {
  wrap_plots(
    c("1-wpi", "2-wpi") |>
      map(function(tp) {
        subset(
          so.inf,
          timepoint == tp &
            celltype %in% c("Monocyte/M\u0278", "Dendritic cell")
        ) |>
          FindMarkers(
            "PINK1 KO", "WT",
            group.by = "genotype", log2fc.threshold = 0.25
          ) |>
          as.data.table(keep.rownames = "gene") |>
          _[, `:=`(
            timepoint = tp
          )][]
      }) |>
      reduce(rbind) |>
      ggplot() +
      aes(
        x = avg_log2FC,
        y = -log(p_val_adj, 10),
        label = fifelse(
          gene %in% c(
            tail(gene[order(abs(avg_log2FC))], 20),
            head(gene[order(abs(p_val_adj))], 20)
          ),
          gene,
          NA
        )
      ) +
      geom_point() +
      geom_text() +
      facet_wrap(vars(timepoint), ncol = 1),
    c("Monocyte", "M\u0278", "Dendritic cell") |>
      map(function(ct) {
        subset(so.inf, celltype.fine == ct) |>
          FindMarkers(
            "PINK1 KO", "WT",
            group.by = "genotype", log2fc.threshold = 0.25
          ) |>
          as.data.table(keep.rownames = "gene") |>
          _[p_val_adj < 0.05, .(
            sign = avg_log2FC > 0,
            gene, celltype = ct
          )]
      }) |>
      reduce(rbind) |>
      _[, .(n = .N), by = .(sign, celltype)] |>
      ggplot() +
      aes(x = 0, y = sign, fill = fifelse(sign, n, -n), label = n) +
      geom_tile() +
      geom_text() +
      facet_wrap(
        vars(
          fct_relevel(celltype, "Monocyte", "M\u0278", "Dendritic cell")
        ),
        ncol = 1
      ) +
      theme_void() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
      guides(fill = "none"),
    ncol = 2,
    widths = c(1, 0.5)
  )
}

fig.4.b <- function(mon.so, mon.pt, plot.embedding = "pca") {
  wrap_plots(
    ... = append(
      FeaturePlot(
        mon.so, c("Ly6c2", "Cx3cr1", "Itgax"),
        reduction = plot.embedding, combine = FALSE
      ),
      list(
        ptp(mon.pt, plot.embedding,
          subset.expr = lineage == "1"
        ) +
          guides(fill = "none") +
          theme_void() +
          theme(strip.text = element_blank()),
        ptp(mon.pt, plot.embedding,
          subset.expr = lineage == "1",
          non.lineage.facet = "genotype",
          show.curves = FALSE, label.by = NULL
        ) +
          theme_void(),
        ptp(mon.pt, plot.embedding,
          subset.expr = lineage == "2"
        ) +
          guides(fill = "none") +
          theme_void() +
          theme(strip.text = element_blank()),
        ptp(mon.pt, plot.embedding,
          subset.expr = lineage == "2",
          non.lineage.facet = "genotype",
          show.curves = FALSE, label.by = NULL
        ) +
          theme_void()
      )
    ),
    design = "
      ABC
      DEE
      FGG
    "
  )
}

fig.4.c <- function(mon.pt) {
  # only looking at infected cells
  mon.stats <- mon.pt[["pseudotime"]][!is.na(pseudotime)][
    mon.pt[["input.meta"]][
      condition == "inf." & celltype.fine == "Monocyte",
      .(genotype, cell, celltype.fine)
    ],
    .(
      mature = pseudotime >= .SD[, mean(pseudotime)],
      genotype, celltype.fine, lineage
    ),
    on = .(cell)
  ][
    ,
    .(n = .N),
    by = .(genotype, celltype.fine, mature, lineage)
  ]
  wrap_plots(
    ggplot(mon.stats) +
      aes(
        x = fct_relevel(genotype, "WT"),
        y = n,
        fill = fct_relevel(fifelse(mature, "mature", "immature"), "immature")
      ) +
      geom_col(position = position_fill(reverse = FALSE)) +
      facet_wrap(
        vars(
          fct_relevel(
            fifelse(lineage == "1", "macrophage-like", "DC-like"),
            "macrophage-like", "DC-like"
          )
        ),
        ncol = 1
      ) +
      labs(y = "proportion", x = NULL, fill = NULL) +
      theme(legend.position = "left"),
    mon.stats[mature == TRUE][
      mon.stats[mature == FALSE],
      .(ratio = n / i.n, genotype, lineage),
      on = .(genotype, lineage)
    ] |>
      ggplot() +
      aes(
        x = fct_relevel(genotype, "WT"),
        y = ratio
      ) +
      geom_col() +
      facet_wrap(
        vars(
          fct_relevel(
            fifelse(lineage == "1", "macrophage-like", "DC-like"),
            "macrophage-like", "DC-like"
          )
        ),
        ncol = 1
      ) +
      labs(y = "ration of mature/immature cells", x = NULL, fill = NULL) +
      guides(fill = "none")
  )
}

fig.4.d <- function(mon.eGO) {
  treeplot(mon.eGO, showCategory = 40)
}

fig.4.e <- function(mon.eGO) {
  cnetplot(mon.eGO, showCategory = 40, max.overlaps = Inf, node_label = "gene")
}

fig.5.a <- function(wt.cc, ko.cc) {
  wrap_plots(
    netAnalysis_signalingRole_scatter(wt.cc),
    netAnalysis_signalingRole_scatter(ko.cc)
  )
}

fig.5.b1 <- function(wt.cc) {
  netVisual_chord_gene(wt.cc,
    sources.use = "Monocyte/M\u0278",
    targets.use = c(
      "B cell", "CD4 Th17", "Monocyte/M\u0278",
      "CD4 Th1", "CD8 Tc1", "CD8 Tn", "Granulocyte",
      "CD4 Th2", "gdT cell", "Dendritic cell", "NK/ILC",
      "dnT cell", "CD4 Treg", "CD4 Tn", "CD4 Trm", "CD8 Trm"
    ),
    slot.name = "netP",
    legend.pos.x = 5, legend.pos.y = 5,
    lab.cex = 2, small.gap = 3,
    big.gap = 20
  )
}

fig.5.b2 <- function(ko.cc) {
  netVisual_chord_gene(ko.cc,
    sources.use = "Monocyte/M\u0278",
    targets.use = c(
      "B cell", "CD4 Th17", "Monocyte/M\u0278",
      "CD4 Th1", "CD8 Tc1", "CD8 Tn", "Granulocyte",
      "CD4 Th2", "gdT cell", "Dendritic cell", "NK/ILC",
      "dnT cell", "CD4 Treg", "CD4 Tn", "CD4 Trm", "CD8 Trm"
    ),
    slot.name = "netP",
    legend.pos.x = 5, legend.pos.y = 5,
    lab.cex = 2, small.gap = 3,
    big.gap = 20
  )
}


# generate figures
if (!interactive()) {
  so <- readRDS("out/so/all.seurat.rds")
  fig <- list()

  fig[["1c"]] <- fig.1.c(so)
  fig[["1d"]] <- fig.1.d(so)
  fig[["1e"]] <- fig.1.e(so)

  fig[["2a"]] <- fig.2.a(so)
  fig[["2b"]] <- fig.2.b(so)

  rm(so)

  fig[["3a"]] <- fig.3.a(
    readRDS("out/pt/cd4.list.rds"),
    readRDS("out/pt/cd8.list.rds")
  )
  fig[["3b"]] <- fig.3.b(
    readRDS("out/so/cd4-inf.seurat.rds"),
    readRDS("out/so/cd8-inf.seurat.rds")
  )
  fig[["3c"]] <- fig.3.c(readRDS("out/eGO/th1.enrichresult.rds"))
  fig[["3d"]] <- fig.3.d(readRDS("out/eGO/tc1.enrichresult.rds"))

  mon.pt <- readRDS("out/pt/mon.list.rds")
  mon.ego <- readRDS("out/eGO/mon.enrichresult.rds")
  fig[["4a"]] <- fig.4.a(readRDS("out/so/all-inf.seurat.rds"))
  fig[["4b"]] <- fig.4.b(readRDS("out/so/mon.seurat.rds"))
  fig[["4c"]] <- fig.4.c(mon.pt)
  fig[["4d"]] <- fig.4.d(mon.ego)
  fig[["4e"]] <- fig.4.e(mon.ego)

  rm(mon.pt, mon.ego)

  wt.cc <- readRDS("out/cc/wt.cellchat.rds")
  ko.cc <- readRDS("out/cc/ko.cellchat.rds")
  fig[["5a"]] <- fig.5.a(wt.cc, ko.cc)
  fig[["5b1"]] <- fig.5.b1(wt.cc)
  fig[["5b2"]] <- fig.5.b2(ko.cc)

  rm(wt.cc, ko.cc)
}

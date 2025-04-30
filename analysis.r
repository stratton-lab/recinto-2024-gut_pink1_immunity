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

rsp <- function(
    so, group.by,
    base.embedding = "pca",
    extra.embeddings = "umap", ...) {
  spt <- slingshot(Embeddings(so, base.embedding), so[[]][[group.by]], ...)

  append(
    setNames(
      list(
        spt,
        group.by,
        slingPseudotime(spt) |>
          as.data.table(keep.rownames = "cell") |>
          melt(
            measure.vars = measure(
              lineage = as.integer,
              pattern = "Lineage([0-9]+)"
            ),
            value.name = "pseudotime"
          ),
        map(c(base.embedding, extra.embeddings), function(emb) {
          `colnames<-`(
            Embeddings(so, emb),
            sprintf("%s_%s", emb, seq_len(ncol(Embeddings(so, emb))))
          )
        }) |>
          purrr::reduce(cbind) |>
          as.data.table(keep.rownames = "cell"),
        as.data.table(so[[]], keep.rownames = "cell"),
        slingCurves(spt, as.df = TRUE) |>
          as.data.table() |>
          `colnames<-`(
            c(
              sprintf("dim_%s", seq_len(ncol(Embeddings(so, base.embedding)))),
              "ord", "lineage"
            )
          ) |>
          _[order(ord)]
      ),
      c(
        "spt", "group.by", "pseudotime",
        "input.embeddings", "input.meta",
        sprintf("curves.%s", base.embedding)
      )
    ),
    map(
      setNames(extra.embeddings, sprintf("curves.%s", extra.embeddings)),
      function(emb) {
        Embeddings(so, emb) |>
          embedCurves(spt, newDimRed = _) |>
          slingCurves(as.df = TRUE) |>
          as.data.table() |>
          `colnames<-`(
            c(
              sprintf("dim_%s", seq_len(ncol(Embeddings(so, emb)))),
              "ord", "lineage"
            )
          ) |>
          _[order(ord)]
      }
    )
  )
}

run_cc <- function(so, group.by, db) {
  cc <- createCellChat(so, group.by = group.by)
  cc@DB <- db
  subsetData(cc) |>
    identifyOverExpressedGenes() |>
    identifyOverExpressedInteractions() |>
    computeCommunProb() |>
    computeCommunProbPathway() |>
    aggregateNet() |>
    netAnalysis_computeCentrality()
}

run_ego <- function(so) {
  FindMarkers(
    so,
    "PINK1 KO",
    "WT",
    group.by = "genotype",
    log2fc.threshold = 0.25
  ) |>
    subset(p_val_adj < 0.05) |>
    rownames() |>
    enrichGO(
      OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP",
      pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.2,
      minGSSize = 2, maxGSSize = 500
    ) |>
    pairwise_termsim(showCategory = 40)
}

sp <- function(so) {
  so |>
    NormalizeData() |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA() |>
    RunUMAP(dims = 1:30)
}

if (TRUE) {
  meta.dt <- fread("data/meta.tsv")
  so <- list(
    "D6/WTuninf" = list(
      "timepoint" = "1-wpi", "genotype" = "WT", "condition" = "uninf."
    ),
    "D13/WTuninf" = list(
      "timepoint" = "2-wpi", "genotype" = "WT", "condition" = "uninf."
    ),
    "D6/WTinf" = list(
      "timepoint" = "1-wpi", "genotype" = "WT", "condition" = "inf."
    ),
    "D13/WTinf" = list(
      "timepoint" = "2-wpi", "genotype" = "WT", "condition" = "inf."
    ),
    "D6/KOuninf" = list(
      "timepoint" = "1-wpi", "genotype" = "PINK1 KO", "condition" = "uninf."
    ),
    "D13/KOuninf" = list(
      "timepoint" = "2-wpi", "genotype" = "PINK1 KO", "condition" = "uninf."
    ),
    "D6/KOinf" = list(
      "timepoint" = "1-wpi", "genotype" = "PINK1 KO", "condition" = "inf."
    ),
    "D13/KOinf" = list(
      "timepoint" = "2-wpi", "genotype" = "PINK1 KO", "condition" = "inf."
    )
  ) |>
    imap(function(info, path) {
      Read10X(sprintf("data/%s/filtered_feature_bc_matrix", path)) |>
        CreateSeuratObject() |>
        RenameCells(add.cell.id = path) |>
        AddMetaData(
          data.frame(
            meta.dt[
              info,
              on = .(condition, genotype, timepoint)
            ][
              , `:=`(cell = sprintf("%s_%s-1", path, cell))
            ][],
            row.names = "cell"
          )
        )
    }) |>
    purrr::reduce(merge) |>
    subset(condition %in% c("inf.", "uninf.")) |> # keep only annotated cells
    JoinLayers() |>
    sp()

  # fix coarse grained annotations using fine grained annotations
  so <- AddMetaData(
    so,
    data.table(
      celltype.fine = c(
        "Mɸ", "Dendritic cell", "Monocyte",
        "CD8 Tn", "CD4 Th2", "CD4 Tn", "CD4 Th2(IL17+)", "CD4 Th1",
        "CD4 Treg", "CD8 Tc1", "CD4 Trm", "CD8 Trm", "CD4 Th17"
      ),
      new.annot = c(
        "Monocyte/Mɸ", "Dendritic cell", "Monocyte/Mɸ",
        "CD8 Tn", "CD4 Th2", "CD4 Tn", "CD4 Th2", "CD4 Th1",
        "CD4 Treg", "CD8 Tc1", "CD4 Trm", "CD8 Trm", "CD4 Th17"
      )
    )[
      as.data.table(so[[c("celltype.fine", "celltype")]], keep.rownames = TRUE),
      .(rn, celltype = fcoalesce(new.annot, celltype)),
      on = .(celltype.fine)
    ] |> data.frame(row.names = "rn")
  )
  # remove CD4/8 prefix from celltype.fine
  so <- AddMetaData(
    so,
    str_remove(so[["celltype.fine"]][, ], "CD(4|8) "),
    "celltype.fine"
  )
}

if (TRUE) {
  so.inf <- subset(so, condition == "inf.")
  cd4.so <- subset(
    so,
    cells = rownames(
      subset(
        so[[]],
        timepoint == "2-wpi" & startsWith(celltype, "CD4")
      )
    )
  ) |> sp()
  cd4.inf.so <- subset(cd4.so, condition == "inf.") |> sp()
  cd8.so <- subset(
    so,
    cells = rownames(
      subset(
        so[[]],
        timepoint == "2-wpi" & startsWith(celltype, "CD8")
      )
    )
  ) |> sp()
  cd8.inf.so <- subset(cd8.so, condition == "inf.") |> sp()
  mon.so <- subset(
    so,
    celltype.fine %in% c("Monocyte", "M\u0278", "Dendritic cell") &
      timepoint == "1-wpi"
  ) |> sp()
  mon.inf.so <- subset(mon.so, condition == "inf.") |> sp()

  cd4.pt <- rsp(cd4.so, "celltype.fine", start.clus = "Tn")
  cd8.pt <- rsp(cd8.so, "celltype.fine", start.clus = "Tn")
  mon.pt <- rsp(mon.so, "celltype.fine", start.clus = "Monocyte")
}

if (TRUE) {
  so[["cd4.pca"]] <- cd4.so[["pca"]]
  so[["cd8.pca"]] <- cd8.so[["pca"]]
  so[["mon.pca"]] <- mon.so[["pca"]]
  so[["cd4.umap"]] <- cd4.so[["umap"]]
  so[["cd8.umap"]] <- cd8.so[["umap"]]
  so[["mon.umap"]] <- mon.so[["umap"]]
  df <- cd4.pt$pseudotime |>
    dcast(cell ~ lineage, value.var = "pseudotime") |>
    data.frame(row.names = "cell")
  colnames(df) <- sprintf("cd4.pseudotime.%s", str_remove(colnames(df), "X"))
  so <- AddMetaData(so, df)
  df <- cd8.pt$pseudotime |>
    dcast(cell ~ lineage, value.var = "pseudotime") |>
    data.frame(row.names = "cell")
  colnames(df) <- sprintf("cd8.pseudotime.%s", str_remove(colnames(df), "X"))
  so <- AddMetaData(so, df)
  df <- mon.pt$pseudotime |>
    dcast(cell ~ lineage, value.var = "pseudotime") |>
    data.frame(row.names = "cell")
  colnames(df) <- sprintf("mon.pseudotime.%s", str_remove(colnames(df), "X"))
  so <- AddMetaData(so, df)

  list(
    "all" = so,
    "all-inf" = so.inf,
    "cd4" = cd4.so,
    "cd4-inf" = cd4.inf.so,
    "cd8" = cd8.so,
    "cd8-inf" = cd8.inf.so,
    "mon" = mon.so,
    "mon-inf" = mon.inf.so,
  ) |> iwalk(function(x, i) {
    saveRDS(sprintf("out/so/%s.seurat.rds", i), x)
  })
  list(
    "cd4" = cd4.pt,
    "cd8" = cd8.pt,
    "mon" = mon.pt
  ) |> iwalk(function(x, i) {
    saveRDS(sprintf("out/pt/%s.list.rds", i), x)
  })
  fwrite(cd4.pt$curves.pca, "out/curves/cd4-curves-pca.csv")
  fwrite(cd4.pt$curves.umap, "out/curves/cd4-curves-umap.csv")
  fwrite(cd8.pt$curves.pca, "out/curves/cd8-curves-pca.csv")
  fwrite(cd8.pt$curves.umap, "out/curves/cd8-curves-umap.csv")
  fwrite(mon.pt$curves.pca, "out/curves/mon-curves-pca.csv")
  fwrite(mon.pt$curves.umap, "out/curves/mon-curves-umap.csv")
}

if (TRUE) {
  th1.eGO <- run_ego(subset(cd4.inf.so, celltype.fine == "Th1"))
  tc1.eGO <- run_ego(subset(cd8.inf.so, celltype.fine == "Tc1"))
  mon.eGO <- run_ego(subset(mon.inf.so, celltype.fine == "Monocyte"))

  list(
    "th1" = th1.eGO,
    "tc1" = tc1.eGO,
    "mon" = mon.eGO
  ) |> iwalk(function(x, i) {
    saveRDS(sprintf("out/eGO/%s.enrichresult.rds", i), x)
  })
}

if (TRUE) {
  wt.cc <- run_cc(
    subset(
      so,
      timepoint == "1-wpi" & condition == "inf." & genotype == "WT"
    ),
    "celltype", CellChatDB.mouse
  )
  ko.cc <- run_cc(
    subset(
      so,
      timepoint == "1-wpi" & condition == "inf." & genotype == "PINK1 KO"
    ),
    "celltype", CellChatDB.mouse
  )

  list(
    "wt" = wt.cc,
    "ko" = ko.cc
  ) |> iwalk(function(x, i) {
    saveRDS(sprintf("out/cc/%s.cellchat.rds", i), x)
  })
}

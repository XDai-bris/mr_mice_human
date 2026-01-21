# traitList <- unique(miceMR$exp)
# all_gwas <- ieugwasr::gwasinfo()
# traitHuman <- all_gwas$trait
# tmp <- cbind(all_gwas$trait, all_gwas$id)
# write.table(tmp, "openGWAS_huamnTraits.csv", row.names = F)
# 
# traitHuman[grepl("fat", traitHuman)]
# tmp[grepl("Whole body fat-free mass", traitHuman),]
# c <- all_gwas[grepl("Total body bone mineral density", traitHuman),]; View(c)
# ---- Libraries ----
source("./FUN.R")
# test OpenGWAS connection
api_status()
Sys.getenv("R_ENVIRON_USER")
ieugwasr::get_opengwas_jwt()
user()

miceMR   <- read.csv("./res_mr_SF_norm.csv")
# phen_map <- read.csv("./micePheno.csv")

res <- mouse_human_mr_onepair(
  exp_mice = "body_length",
  out_mice = "uCT_BMD",
  human_exp_id = "ukb-a-389",
  human_outcome_id = "ebi-a-GCST005348",
  miceMR = miceMR,
  human_exp_label = "Standing height",
  human_out_label = "Total body bone mineral density",
  make_plot = TRUE
)

print(res)         # neat console summary
res$summary        # tidy one-row tibble for downstream
if (!is.null(res$plot)) print(res$plot)  # forest plot

# run all pair of mice-human MR

# micePheno <- read.csv("./micePheno.csv", header = T)

# unique_rows <- micePheno[!duplicated(micePheno$openGWASID), ]
# write.csv(unique_rows, file = "inputOpenGWAS.csv", row.names = F)
rows <- read.csv("./inputOpenGWAS.csv", header = T)
n <- nrow(rows)
# columns in your unique_rows
trait_col <- "Phenotype"      # mouse traits
label_col <- "HumanTrait"     # human labels
id_col <- "openGWASID"        # human GWAS IDs

# helper to make safe list names
sanitize <- function(x) {
  x <- gsub("\\s+", "_", x)
  gsub("[^A-Za-z0-9_\\-]+", "", x)
}

# all ordered pairs (i != j)
pairs <- subset(expand.grid(i = seq_len(n), j = seq_len(n)), i != j)

# preallocate results
results <- vector("list", nrow(pairs))
names(results) <- sprintf(
  "exp:%s_out:%s",
  sanitize(rows[[trait_col]][pairs$i]),
  sanitize(rows[[trait_col]][pairs$j])
)

errors <- list()

for (r in seq_len(nrow(pairs))) {
  i <- pairs$i[r]; j <- pairs$j[r]
  
  exp_mice  <- rows[[trait_col]][i]
  out_mice  <- rows[[trait_col]][j]
  human_exp_id  <- rows[[id_col]][i]
  human_out_id  <- rows[[id_col]][j]
  
  human_exp_label <- rows[[label_col]][i]
  human_out_label <- rows[[label_col]][j]
  
  nm <- names(results)[r]
  message("Running: ", nm)
  
  res <- try(
    mouse_human_mr_onepair(
      exp_mice = exp_mice,
      out_mice = out_mice,
      human_exp_id = human_exp_id,
      human_outcome_id = human_out_id,
      miceMR = miceMR,
      human_exp_label = human_exp_label,
      human_out_label = human_out_label,
      make_plot = TRUE
    ),
    silent = TRUE
  )
  
  if (inherits(res, "try-error")) {
    errors[[nm]] <- as.character(res)
    results[[r]] <- NULL
  } else {
    results[[r]] <- res
  }
}

# combined summary
summary_tbl <- bind_rows(
  lapply(results[!vapply(results, is.null, logical(1))], function(x) {
    x$summary %>%
      mutate(run_name = paste0("exp:", x$mapping$mouse_exp,
                               "_out:", x$mapping$mouse_out))
  })
)

# check counts
cat("Total possible runs:", nrow(pairs), "\n")
cat("Successful runs:", sum(!vapply(results, is.null, logical(1))), "\n")
cat("Failed runs:", length(errors), "\n")

# optionally save
# saveRDS(results, "mouse_human_mr_grid_results.rds")
write.csv(summary_tbl, "mice_human_mr.csv", row.names = FALSE)

summary_tbl_flt <- summary_tbl[which(summary_tbl$human_nsnp > 10), ]
write.csv(summary_tbl_flt, "mice_human_mr_sum_flt.csv", row.names = FALSE)


# Historical prototype: original exploratory SNP/window scan retained for
# reference only. It predates the staged workflow and is not the recommended
# entry point for routine runs or parameter tuning.
# Inputs/outputs may not match config/config.yml or the numbered scripts.

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(Biostrings)
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

setwd("C:/Users/jc213439/Dropbox/Emmanuel/Whitespot/hvr_scan")

# -----------------------------
# Hard-code inputs
# -----------------------------
fasta_path <- "../whole_genome_alignments/41seqscollagen_Edited.fasta"
tree_path  <- "./data/41SEQSUPDATEDEDITEDMCC"   # or .nex/.nexus


outdir <- "wssv_all_out"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read alignment
# -----------------------------
dna <- readDNAStringSet(fasta_path)
if (length(dna) < 2) stop("Alignment FASTA contains <2 sequences.")
aln_names <- names(dna)
if (anyDuplicated(aln_names)) stop("Duplicate sequence names in FASTA.")

seq_lengths <- width(dna)
if (length(unique(seq_lengths)) != 1) {
  stop("Sequences have different lengths. Ensure this is an alignment (all same length).")
}
L <- unique(seq_lengths)

# -----------------------------
# Read tree
# -----------------------------
read_tree_any <- function(path) {
  tr <- tryCatch(read.tree(path), error = function(e) NULL)
  if (!is.null(tr)) return(tr)
  tr2 <- tryCatch(read.nexus(path), error = function(e) NULL)
  if (!is.null(tr2)) return(tr2)
  stop("Failed to read tree as Newick or Nexus: ", path)
}
tree <- read_tree_any(tree_path)

if (is.null(tree$tip.label) || length(tree$tip.label) < 2) stop("Tree has <2 tips.")
if (!is.binary(tree)) {
  message("Tree is not strictly binary; consider resolving polytomies.")
}

# -----------------------------
# Match labels and prune/reorder
# -----------------------------
common <- intersect(tree$tip.label, aln_names)

cat("Alignment: n =", length(dna), " L =", L, "\n")
cat("Tree:      tips =", length(tree$tip.label), " internal =", tree$Nnode, "\n")
cat("Shared labels between tree and alignment:", length(common), "\n")

if (length(common) < 2) stop("Too few shared labels between tree tips and alignment names.")

tree2 <- drop.tip(tree, setdiff(tree$tip.label, common))
dna2  <- dna[tree2$tip.label]  # reorder to match tree tips
stopifnot(identical(names(dna2), tree2$tip.label))

tree <- tree2
dna <- dna2


# -----------------------------
# Extract recombinant-free tails
# -----------------------------
extract_and_merge_windows <- function(dna, windows) {
  pieces <- lapply(windows, \(r) subseq(dna, start = r[1], end = r[2]))
  dna_slice <- Reduce(function(a,b) xscat(a,b), pieces)
  names(dna_slice) <- names(dna)
  return(dna_slice)
}

tail_sites <- list(
  right = c(start=297825,end=L),
  left = c(start=1,end=16056)
)


tails <- extract_and_merge_windows(dna, tail_sites)

cat("Extracted tails: n =", length(tails), " L =", width(tails)[1], "\n")


# -----------------------------
# 1) Integer encoding (A,C,G,T -> 1..4; everything else -> 0)
# -----------------------------
encode_alignment_int <- function(dna_set) {
  stopifnot(inherits(dna_set, "DNAStringSet"))
  N <- length(dna_set)
  L <- unique(width(dna_set))
  if (length(L) != 1) stop("encode_alignment_int: sequences have unequal lengths")
  
  # ASCII lookup table: raw byte -> integer code
  lut <- integer(256)
  lut[as.integer(charToRaw("A")) + 1L] <- 1L
  lut[as.integer(charToRaw("C")) + 1L] <- 2L
  lut[as.integer(charToRaw("G")) + 1L] <- 3L
  lut[as.integer(charToRaw("T")) + 1L] <- 4L
  lut[as.integer(charToRaw("a")) + 1L] <- 1L
  lut[as.integer(charToRaw("c")) + 1L] <- 2L
  lut[as.integer(charToRaw("g")) + 1L] <- 3L
  lut[as.integer(charToRaw("t")) + 1L] <- 4L
  
  aln <- matrix(0L, nrow = N, ncol = L)
  for (i in seq_len(N)) {
    r <- charToRaw(as.character(dna_set[[i]]))
    aln[i, ] <- lut[as.integer(r) + 1L]
  }
  rownames(aln) <- names(dna_set)
  aln
}

# aln_int <- encode_alignment_int(tails)
aln_int <- encode_alignment_int(dna)
Ntip <- length(tree$tip.label)
L <- ncol(aln_int)

cat("Encoded alignment int matrix:", nrow(aln_int), "x", ncol(aln_int), "\n")
cat("Missing/ambiguous rate:", mean(aln_int == 0), "\n")


# -----------------------------
# 2) Per-node split masks (left side defined by first child clade)
# -----------------------------
get_children <- function(tree, node) {
  tree$edge[tree$edge[, 1] == node, 2]
}

get_tip_desc <- function(tree, node) {
  if (node <= length(tree$tip.label)) return(node)
  phangorn::Descendants(tree, node, type = "tips")[[1]]
}

precompute_node_clades <- function(tree, min_clade_size = 2L, max_clade_frac = 0.95) {
  Ntip <- length(tree$tip.label)
  internal_nodes <- (Ntip + 1):(Ntip + tree$Nnode)
  
  root <- Ntip + 1L
  candidate_nodes <- setdiff(internal_nodes, root)
  
  clade_mask <- matrix(FALSE, nrow = length(candidate_nodes), ncol = Ntip)
  node_id <- integer(length(candidate_nodes))
  clade_tips <- vector("list", length(candidate_nodes))
  
  keep <- logical(length(candidate_nodes))
  
  for (i in seq_along(candidate_nodes)) {
    node <- candidate_nodes[i]
    tips <- get_tip_desc(tree, node)
    
    n_clade <- length(tips)
    
    if (
      n_clade >= min_clade_size &&
      n_clade <= floor(max_clade_frac * Ntip)
    ) {
      clade_mask[i, tips] <- TRUE
      clade_tips[[i]] <- tips
      node_id[i] <- node
      keep[i] <- TRUE
    }
  }
  
  list(
    node_id = node_id[keep],
    clade_tips = clade_tips[keep],
    clade_mask = clade_mask[keep, , drop = FALSE]
  )
}

clades <- precompute_node_clades(tree)
target_mask <- clades$clade_mask


# -----------------------------
# Node info
# -----------------------------
node_depth_edges <- ape::node.depth(tree2)  # length Ntip + Nnode
depth_internal <- node_depth_edges[splits$node_id]

n_left  <- rowSums(splits$left_mask)
n_right <- Ntip - n_left

# balance in [0,1], max at perfectly balanced split
balance <- (2 * pmin(n_left, n_right)) / Ntip

# depth weight: near-root higher (depth=0 => 1, depth=1 => 1/2, ...)
depth_w <- 1 / (depth_internal + 1)

# combined weight (tweakable)
node_weight <- depth_w * balance

node_info <- data.frame(
  node = splits$node_id,
  depth = depth_internal,
  n_left = n_left,
  n_right = n_right,
  balance = balance,
  weight = node_weight
)

cat("Internal nodes:", nrow(node_info), "\n")
cat("Weight summary:\n")
print(summary(node_info$weight))

# outdir <- "wssv_tails_out"  # keep consistent with Step 1
saveRDS(
  list(
    tree = tree2,
    aln_int = aln_int,
    splits = splits,
    node_info = node_info
  ),
  file = file.path(outdir, "precomputed.rds")
)



# -----------------------------
# Read in preprocessed outputs
# -----------------------------
pre <- readRDS(file.path(outdir, "precomputed.rds"))

aln_int   <- pre$aln_int            # Ntip x L, values 0..4
target_mask <- pre$splits$target_mask   # Nnode x Ntip (logical)
node_info <- pre$node_info

Ntip  <- nrow(aln_int)
L     <- ncol(aln_int)
Nnode <- nrow(target_mask)

node_weight <- node_info$weight


# -----------------------------
# Score site / allele predictability
# -----------------------------
min_total_obs  <- 30L  # require at least this many observed tips at a node to consider it for scoring
min_side_obs   <- 2L   # require at least this many observed tips on each side of the node to consider it for scoring
min_site_maf   <- 2L   # require at least this many observed tips with the minor allele at this site to consider it for scoring

# -----------------------------
# Same best-rule function, but ALSO return obs_left/right/total_obs
# -----------------------------
best_rule_for_site <- function(x, target_mask, min_total_obs, min_side_obs, min_site_maf) {
  obs <- (x != 0L)
  tab_all <- tabulate(x[obs], nbins = 4L)
  
  obs_in  <- as.vector(target_mask %*% as.integer(obs))
  obs_out <- as.vector((!target_mask) %*% as.integer(obs))
  total_obs <- obs_in + obs_out
  
  ok_node <- (total_obs >= min_total_obs) &
    (obs_in >= min_side_obs) &
    (obs_out >= min_side_obs)
  
  best_acc    <- rep(NA_real_, nrow(target_mask))
  best_allele <- rep(NA_integer_, nrow(target_mask))
  best_dir    <- rep(NA_integer_, nrow(target_mask))
  
  best_acc[ok_node] <- -Inf
  
  for (a in 1:4) {
    if (tab_all[a] < min_site_maf) next
    
    xa <- (x == a)
    
    c_in  <- as.vector(target_mask %*% as.integer(xa))
    c_out <- as.vector((!target_mask) %*% as.integer(xa))
    
    acc_marks_clade <- (c_in + (obs_out - c_out)) / total_obs
    acc_marks_outside <- (c_out + (obs_in - c_in)) / total_obs
    
    acc_a <- pmax(acc_marks_clade, acc_marks_outside)
    dir_a <- ifelse(acc_marks_clade >= acc_marks_outside, 1L, 2L)
    
    improve <- ok_node & (acc_a > best_acc)
    
    if (any(improve)) {
      best_acc[improve] <- acc_a[improve]
      best_allele[improve] <- a
      best_dir[improve] <- dir_a[improve]
    }
  }
  
  bad <- ok_node & !is.finite(best_acc)
  best_acc[bad] <- NA_real_
  
  list(
    best_acc = best_acc,
    best_allele = best_allele,
    best_dir = best_dir,  # 1 = allele marks clade, 2 = allele marks outside
    obs_in = obs_in,
    obs_out = obs_out,
    total_obs = total_obs
  )
}




# -----------------------------
# Scan all and score all sites
# -----------------------------

# Filter for polymorphic sites
is_snp_col_int <- function(x) {
  x_obs <- x[x != 0L]          # keep only observed A/C/G/T
  length(unique(x_obs)) > 1L   # polymorphic among observed
}

# count a node as "helped" if site improves accuracy by at least this much
min_gain       <- 1/Ntip  # = 1/40; one additional correct tip if total_obs ~ 40
chunk_size     <- 5000L

scores_raw <- matrix(-Inf, nrow=Nnode, ncol=L)  # global gain score for each site
scores_norm <- matrix(-Inf, nrow=Nnode, ncol=L)  # global gain score for each site


cat("Scanning", L, "sites with gain scoring...\n")

for (start in seq(1L, L, by = chunk_size)) {
  end <- min(L, start + chunk_size - 1L)
  
  for (s in start:end) {
    x <- aln_int[, s]
    if (!is_snp_col_int(x)) next  # skip non-polymorphic sites
    br <- best_rule_for_site(
      x = x,
      target_mask = target_mask,
      min_total_obs = min_total_obs,
      min_side_obs = min_side_obs,
      min_site_maf = min_site_maf
    )
    
    ok <- !is.na(br$best_acc)
    if (!any(ok)) next
    
    baseline <- pmax(br$obs_in, br$obs_out) / br$total_obs
    gain <- br$best_acc - baseline
    
    max_gain <- 1 - baseline
    gain_norm <- gain / max_gain
    gain_norm[!is.finite(gain_norm)] <- NA_real_
    
    min_gain_norm <- 0.5
    
    helped <- ok &
      (gain >= 1 / br$total_obs) &
      (gain_norm >= min_gain_norm) &
      is.finite(gain_norm)
    if (!any(helped)) next
    
    gain[!helped] <- -Inf
    gain_norm[!helped] <- -Inf
    
    scores_raw[, s]  <- gain
    scores_norm[, s] <- gain_norm
    
  }
  
  if ((start == 1L) || (end %% (10L * chunk_size) == 0L) || (end == L)) {
    cat("  processed sites", start, "to", end, "\n")
  }
}

saveRDS(
  list(
    params = list(
      min_total_obs = min_total_obs,
      min_side_obs = min_side_obs,
      min_site_maf = min_site_maf,
      min_gain = min_gain,
      chunk_size = chunk_size
    ),
    scores_raw = scores_raw,
    scores_norm = scores_norm
  ),
  file = file.path(outdir, "all_site_scores.rds")
)


# ------------------------------------------------------------
# 0) Normalize: treat -Inf as NA for analysis
# ------------------------------------------------------------
scores2 <- scores_raw
scores2[!is.finite(scores2)] <- NA_real_

finite_n <- sum(!is.na(scores2))
cat("Finite entries:", finite_n, "out of", Nnode * L, "\n")

scores3 <- scores_norm
scores3[!is.finite(scores3)] <- NA_real_

finite_n <- sum(!is.na(scores3))
cat("Finite entries:", finite_n, "out of", Nnode * L, "\n")

# ------------------------------------------------------------
# 1) Convert to a sparse long table of only finite entries
#    (keeps memory sane and makes summaries easy)
# ------------------------------------------------------------
idx <- which(!is.na(scores3), arr.ind = TRUE)
dt <- data.table(
  node_i = idx[, 1],
  site   = idx[, 2],
  gain   = scores3[idx]
)

# Attach node metadata
dt[, `:=`(
  node   = node_info$node[node_i],
  depth  = node_info$depth[node_i],
  n_left = node_info$n_left[node_i],
  n_right= node_info$n_right[node_i],
  balance= node_info$balance[node_i],
  weight = node_info$weight[node_i]
)]

setorder(dt, node, -gain)

cat("Long table rows (finite gains):", nrow(dt), "\n")

dt %>%
  filter(node %in% c(42, 43, 44, 65, 71, 72, 80)) %>%
  ggplot(aes(x=site, y=gain, col=as.factor(node))) +
  geom_point(alpha=0.5) +
  # geom_jitter() +
  labs(title="Gain scores for node-site pairs",
       x="Site index",
       y="Gain score") +
  theme_minimal()


# ------------------------------------------------------------
# 3) Per-site summaries
#    - how many nodes helped per site
#    - max gain at site
#    - weighted sum of gains across nodes (uses node_info$weight)
# ------------------------------------------------------------
site_summary <- dt[, .(
  nodes_helped = .N,
  gain_max     = max(gain),
  gain_mean    = mean(gain),
  gain_sum     = sum(gain),
  wg_sum       = sum(weight * gain),
  wg_mean      = sum(weight * gain) / sum(weight)
), by = site]

# Rank sites by different criteria
setorder(site_summary, -wg_sum)

cat("Sites with any finite gain:", nrow(site_summary), "\n")

# Top sites table
top_sites <- site_summary[1:min(50, .N)]
print(top_sites)






# ------------------------------------------------------------
# 4) Per-node summaries
#    - how many sites help each node
#    - best site for each node
# ------------------------------------------------------------
node_summary <- dt[, .(
  sites_helpful = .N,
  gain_max      = max(gain),
  gain_mean     = mean(gain),
  best_site     = site[which.max(gain)]
), by = node_i]

node_summary[, `:=`(
  node   = node_info$node[node_i],
  depth  = node_info$depth[node_i],
  n_left = node_info$n_left[node_i],
  n_right= node_info$n_right[node_i],
  balance= node_info$balance[node_i],
  weight = node_info$weight[node_i]
)]

setorder(node_summary, -gain_max)

cat("Nodes with any finite gain:", nrow(node_summary), "\n")
print(node_summary[1:min(30, .N)])




# ------------------------------------------------------------
# 5) Extract “top K” sites per node and “top K” nodes per site
# ------------------------------------------------------------
topK_per_node <- function(dt, K = 10L) {
  dt[, head(.SD, K), by = node_i]
}
topK_per_site <- function(dt, K = 10L) {
  dt[, head(.SD, K), by = site]
}

dt_top_node <- topK_per_node(dt, K = 10L)
dt_top_site <- topK_per_site(dt, K = 10L)


# Just look at the tails
tail_sites <- c(1:16055, 297825:L)
setorder(dt, -gain)
dt_tails <- dt[site %in% tail_sites, ]

dt_top_node <- topK_per_node(dt_tails, K = 10L)
dt_top_site <- topK_per_site(dt_tails, K = 10L)



# ------------------------------------------------------------
# 6) Build PCR candidates
# ------------------------------------------------------------
infer_ml_tree <- function(aln, model = "GTR", k = 4, optNni = TRUE, rearrangement = "stochastic") {
  # aln: character matrix rows=seqs cols=sites; must be A/C/G/T only for this routine.
  if (ncol(aln) < 10) stop("Too few sites after filtering to infer an ML tree.")
  
  phy <- phyDat(aln, type = "DNA")  # requires A/C/G/T
  dm <- dist.ml(phy)
  nj_tree <- NJ(dm)
  fit <- pml(nj_tree, data = phy)
  
  # Topology-focused optimisation
  fit <- optim.pml(
    fit,
    model = model,
    optInv = FALSE,
    optGamma = TRUE,
    k = k,
    optEdge = TRUE,
    optNni = optNni,
    rearrangement = rearrangement
  )
  
  fit$tree
}


# candidate_window <- extract_and_merge_windows(dna, list(c(307919, 308491), # Node 41
#                                                         c(8707, 9893), # Nodes 43, 56
#                                                         c(300765, 302742) # Node 53
#                                                         # c(3952, 4893),  # Nodes 63, 70
#                                                         # c(9893, 10051) # Node 54
#                                                         # c(13186, 13597),
#                                                         # c(6518, 7228),
#                                                         # c(297825, 299571),
# )
# )

candidate_window <- extract_and_merge_windows(dna, list(c(308012, 308972),
                                                        c(19569, 19616),
                                                        c(22650, 23979)
                                                        # c(26417, 26446) # Node 53
                                                        # c(3952, 4893),  # Nodes 63, 70
                                                        # c(9893, 10051) # Node 54
                                                        # c(13186, 13597),
                                                        # c(6518, 7228),
                                                        # c(297825, 299571),
)
)

# dna is a DNAStringSet with aligned sequences
stopifnot(length(unique(width(candidate_window))) == 1)

# Convert to matrix: rows = sequences, cols = sites
candidate_aln_mat <- as.matrix(candidate_window)

# Function to test if a column is polymorphic (A/C/G/T only)
is_snp_col <- function(col) {
  nuc <- col[col %in% c("A","C","G","T")]
  length(unique(nuc)) > 1
}

# Identify SNP columns
snp_cols <- apply(candidate_aln_mat, 2, is_snp_col)

cat("Number of SNP sites:", sum(snp_cols), "\n")


candidate_tree = infer_ml_tree(as.matrix(candidate_window), model = "GTR", k = 4, optNni = TRUE, rearrangement = "stochastic")
plot(candidate_tree)

plot(tree)
nodelabels()

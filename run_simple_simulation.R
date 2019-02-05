library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(multidplyr)

# Simulation function

sim_microbiome <- function(num_subj, 
                           tot_num_reads, 
                           prior.y, 
                           prior.x, 
                           seed_num, 
                           theta = 0.5,
                           sigma_L = 0.77) {
  
  # how many bugs present is not decided by researcher
  num_bug <- nrow(prior.y)
  
  # 
  set.seed(seed_num)
  
  #
  microbe_ids <- paste0("OTU", seq(from = 1, to = num_bug, by = 1))
  subj_ids <- paste0("Subj", seq(from = 1, to = num_subj, by = 1))
  
  # For sample i, the value of covariate (diet) was
  # drawn from a emprical distribution (public data, Z scores)
  x <- sample( prior.x$x, size = num_subj, replace = TRUE )
  names(x) <- subj_ids
  
  # eta[j] = 1 if the covariate influences the abundance of the jth taxa
  # otherwise, eta[j] = 0
  eta <- rbernoulli( num_bug, p = theta )
  names(eta) <- microbe_ids
  
  # the regression parameter b[j] captures 
  # the effect of the covariate on the abundance for taxon j.
  # We used t7 distribution with scale 2.5
  b <- 0 * (1 - eta) + eta * rt.scaled( num_bug, df = 7, mean = 0, sd = 2.5 )
  
  # ilogit(a[i,j]) corresponds to the 
  # baseline abundance for the taxon j in sample i.
  # For sample i, ilogit(a[i,]) is drawn from a 
  # Dirichlet distribution
  ilogit_a <- rdiric( num_subj, shape = prior.y$alpha0 )
  a <- logit( ilogit_a )
  row.names(a) <- subj_ids
  colnames(a) <- microbe_ids
  
  # mu[i,j] represents the abundance of taxon j in sample i
  # logistic regression
  mu <- logit( a + x %*% t(b), inverse = TRUE )
  row.names(mu) <- subj_ids
  colnames(mu) <- microbe_ids
  
  # library sizes were randomly drawn from 
  # a lognormal distribution.
  # This step mimicked varying total reads per sample.
  N <- rlnorm( num_subj, meanlog = log(tot_num_reads), sdlog = sigma_L )
  N <- ceiling(N)
  names(N) <- subj_ids
  
  # For sample i, the observed OTU counts were 
  # randomly generated from a multinomial distribution
  y <- matrix( nrow = num_bug, ncol = num_subj )
  row.names(y) <- microbe_ids
  colnames(y) <- subj_ids
  
  for ( i in 1:num_subj ) {
    y[,i] <- rmultinom( 1, size = N[i], prob = mu[i,] )
  }
  
  output <- list()
  output$y <- y
  output$x <- x
  output$eta <- eta
  return(output)
}

# Analyze functions

## Preprocess otu table
preprocess_otu <- function(microbiome) {
  
  Y <- microbiome$y
  
  # remove OTUs with fewer than 10 reads 
  # as well as OTUs which were present in fewer than 1% of samples
  Y_new1 <- Y[rowSums(Y) > 10,]
  Y_new2 <- Y_new1[rowMeans(Y_new1 != 0) > 0.01,]
  return(Y_new2)
}

# Analyze simulation data

## Zero-inflated Gaussian (metagenomeSeq)

run_metagenomeSeq <- function(microbiome) {
  
  Y <- preprocess_otu(microbiome)
  
  my_obj <- newMRexperiment(counts = Y)
  
  my_obj <- cumNorm(my_obj, p = cumNormStatFast(my_obj))
  x <- microbiome$x
  my_design <- model.matrix(~x)
  my_fit <- fitZig(my_obj, mod = my_design)
  
  my_res <- tidy(my_fit$eb) %>% 
    filter(term == "x") %>% 
    select(gene, estimate, p.value)
  my_res$fdr <- p.adjust(my_res$p.value, method = "BH")
  
  return(my_res)
}

# Poisson regression (shotgunFunctionalizeR)
nest_xNy <- function(microbiome) {
  y_df <- microbiome$y %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    dplyr::rename(gene = rowname) %>% 
    tbl_df() %>% 
    gather(who, y, -gene)
  
  N_df <- y_df %>% 
    group_by(who) %>% 
    summarize(N = sum(y))
  
  x_df <- tibble(who = names(microbiome$x), x = microbiome$x)
  
  y_df %>%  
    left_join(N_df, by = "who") %>% 
    left_join(x_df, by = "who") %>% 
    group_by(gene) %>% 
    nest()
}

add_pois <- function(df) {
  glm(y ~ x, offset = log(N), data = df, family = "poisson")
}

run_shotgunFunctionalizeR <- function(microbiome) {
  
  microbiome$y <- preprocess_otu(microbiome)
  
  y_nest <- nest_xNy(microbiome)
  
  my_res <- y_nest %>% 
    mutate(pois_reg = map(data, add_pois), 
           stats = map(pois_reg, tidy)) %>% 
    dplyr::select(gene, stats) %>% 
    unnest(stats) %>% 
    filter(term == "x") %>% 
    select(gene, estimate, p.value)
  
  my_res$fdr <- p.adjust(my_res$p.value, method = "BH")
  return(my_res)
}

# Voom + limma
run_voomLimma <- function(microbiome) {
  x <- microbiome$x
  X <- model.matrix(~x)
  
  # OTU count is the response variable
  Y <- preprocess_otu(microbiome)
  
  ## voom transformation
  voomY <- limma::voom(Y, design = X)
  
  my_fit <- lmFit(voomY$E, design = voomY$design, weights = voomY$weights)
  my_fit <- eBayes(my_fit)
  my_res <- tidy(my_fit) %>% 
    select(gene, estimate, p.value)
  
  my_res$fdr <- p.adjust(my_res$p.value, method = "BH")
  return(my_res)
  
}

## edgeR

run_edger <- function(microbiome) {
  
  x <- microbiome$x
  X <- model.matrix(~x)
  
  # OTU count is the response variable
  Y <- preprocess_otu(microbiome)
  
  my_dispersions <- estimateDisp(Y, design = X)
  my_fit <- glmFit(Y, 
                   design = X, 
                   dispersion = my_dispersions$tagwise.dispersion, 
                   offset = log(colSums(Y)))
  my_test <- glmLRT(my_fit, coef = "x")
  
  # tidy the outputs
  my_res <- my_test$table %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    rename(gene = rowname, 
           p.value = PValue, 
           estimate = logFC) %>% 
    tbl_df() %>% 
    select(gene, estimate, p.value)
  
  my_res$fdr <- p.adjust(my_res$p.value, method = "BH")
  return(my_res)
}

# Run simulations

run_simulations <- function(num_subj, 
                            tot_num_reads, 
                            prior.y, 
                            prior.x, 
                            seed_num, 
                            theta = 0.5, 
                            sigma_L = 0.77) {
  # step 1: simulate counts
  microbiome <- sim_microbiome(num_subj, 
                               tot_num_reads, 
                               prior.y, 
                               prior.x, 
                               seed_num, 
                               theta, 
                               sigma_L)
  
  # step 2.1 run metagenomeSeq
  metagenomeSeq.out <- run_metagenomeSeq(microbiome)
  
  # step 2.2 run shotgunFunctionalizeR
  shotgunFunctionalizeR.out <- run_shotgunFunctionalizeR(microbiome)
  
  # step 2.3 run Voom + limma
  limma.out <- run_voomLimma(microbiome)
  
  # step 2.4 run edgeR
  edger.out <- run_edger(microbiome)
  
  output <- list()
  output$skewness <- moments::skewness(microbiome$x)
  output$truth <- microbiome$eta
  output$metagenomeSeq <- metagenomeSeq.out
  output$shotgunFunctionalizeR <- shotgunFunctionalizeR.out
  output$limma <- limma.out
  output$edgeR <- edger.out
  
  return(output)
}

## simulation settings
num_sim = 214             # number of simulations
num_cores = 6             # number of cpus
num_subj <- 1000           # number of subjects
sequencing_depth <- 50000 # sequencing depth (pre-decided)

my_seeds <- sample(seq(from = 1, to = 1e5), size = num_sim, replace = FALSE)
dat.sim <- tibble(simulation = 1:num_sim, 
                  seed_num = my_seeds)

# In every simulation, x is drawn from a empirical distribution
## load prior data (y) estimated from Human Microbiome Project (stool samples)
prior.y <- read_rds("HMP_stool_DM_prior.rds")

## load prior data (x) from Nature publication
prior.x <- read_rds("Wu2011_FFQ.rds")

prior.x.nest <- prior.x %>% 
  gather(what.x, x, -who) %>% 
  group_by(what.x) %>% 
  nest()

#HMP.sim <- sim_microbiome(num_subj = 319, tot_num_reads = 4500, prior.y = prior.y, prior.x = prior.x.nest$data[[1]], seed_num = 1)

#write_rds(HMP.sim, "simulated_HMP_V35_stool.rds")

dat.sim$what.x <- prior.x.nest$what.x
dat.sim$data.x <- prior.x.nest$data

xiang_cluster <- create_cluster(cores = num_cores) %>% 
  cluster_library("purrr") %>% 
  cluster_library("tibble") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("tidyr") %>% 
  cluster_library("biobroom") %>% 
  cluster_library("VGAM") %>% 
  cluster_library("metRology") %>% 
  cluster_library("edgeR") %>% 
  cluster_library("limma") %>% 
  cluster_library("metagenomeSeq") %>% 
  cluster_library("moments") %>% 
  cluster_assign_value("sim_microbiome", sim_microbiome) %>% 
  cluster_assign_value("preprocess_otu", preprocess_otu) %>% 
  cluster_assign_value("run_metagenomeSeq", run_metagenomeSeq) %>% 
  cluster_assign_value("nest_xNy", nest_xNy) %>% 
  cluster_assign_value("add_pois", add_pois) %>% 
  cluster_assign_value("run_shotgunFunctionalizeR", run_shotgunFunctionalizeR) %>% 
  cluster_assign_value("run_voomLimma", run_voomLimma) %>% 
  cluster_assign_value("run_edger", run_edger) %>% 
  cluster_assign_value("run_simulations", run_simulations) %>% 
  cluster_assign_value("num_subj", num_subj) %>% 
  cluster_assign_value("sequencing_depth", sequencing_depth) %>%
  cluster_assign_value("prior.y", prior.y)

dat.sim <- dat.sim %>% 
  partition(simulation, cluster = xiang_cluster) %>% 
  mutate(result = map2(seed_num, data.x, ~run_simulations(seed_num = .x, 
                                                          prior.x = .y, 
                                                          num_subj = num_subj, 
                                                          tot_num_reads = sequencing_depth, 
                                                          prior.y = prior.y))) %>% 
  collect() %>% 
  as_tibble()

write_rds(dat.sim, "simulation_Yhmp_Xffq_N1K_S50K.rds")



# Every function ever ####
# General ####
# Plotting palette
cols <- c('#fff65c', '#2bc6ff', '#ba9863', '#52d9c0', '#3c4278',
          '#4c83e0', '#cc1fb2', '#00ffbf', '#1a7878', '#c47900',
          '#0025ba', '#ebb217', '#d4d4ff', '#009e45', '#fbff00',
          '#7b42ff', '#ff0000','#0000ff', '#ffa1a1', '#1e0078',
          '#780000', '#000000')

# LiDAR point cloud processing ####
# See use in: "study.crown.R", "random.crown.R", "clean.R", and "tree.map.R"
random.plot <- function(min.x, max.x, min.y, max.y) {
  # Accepting the min and max dimensions of a laz file as input, a number of
  # randomly sampled xy coordinates are sampled from that file and returned as output. The number of sampled xy pairs
  # is proportional to the area of the file.
  
  # I realize that it's not necessary to create a fraction
  # "file.area/cat.area". This will be adjusted in future applications
  
  # I'm also not a fan of randomly sampling xy pairs @ 1 meter resolution, so this
  # will be updated in future applications. I'd rather the sample resolution match the
  # resolution of the points, but it shouldn't affect anything.
  
  # If you want a certain number of centers to be sampled, you should
  # change the hard-coded "8500" to something else.
  
  cat.x.dim <- 571803.7 - 558288.2
  cat.y.dim <- 4904443 - 4894461
  cat.area <- cat.x.dim * cat.y.dim
  
  x.dim <- max.x - min.x
  y.dim <- max.y - min.y
  file.area <- x.dim * y.dim
  
  centers <- file.area / cat.area * 8500 %>% round() 
  
  x.center <- sample(round(min.x):round(max.x), centers)
  y.center <- sample(round(min.y):round(max.y), centers)
  
  tibble(x.center = x.center,
         y.center = y.center)
}
tree.samp <- function(norm.list, cores){
  # This function accepts a list of laz file names associated with ground-normalized 50 m radius-clips
  # (i.e. all ground points have elevations adjusted to 0, while non-ground points above
  # are adjusted accordingly), applies a tree-detecting algorithm to each clip, randomly chooses a single
  # tree found in the clip, and outputs a data frame containing locations and the maximum height from which
  # randomly selected height can be selected.
  
  require(foreach)
  require(doMC)
  registerDoMC(cores = cores)
  
  set.seed(666)
  samp.trees <- foreach(i = 1:length(norm.list), .combine = rbind, .errorhandling = 'remove') %dopar% {
    
    las <- readLAS(norm.list[[i]])
    
    trees <- find_trees(las,
                        algorithm = lmf(ws = 5, hmin = 35, shape = 'circular'),
                        uniqueness = 'bitmerge')
    
    trees %<>% .@coords %>% data.frame() %>% cbind(trees@data) %>%
      distinct(treeID, .keep_all = T)
    
    ids <- trees$treeID
    
    samp.ids <- sample(ids, size = 1)
    trees %<>% filter(treeID %in% samp.ids)
    
    trees$max.samp <- trees$Z - 5
    trees$max.samp[trees$max.samp > 55] <- 55
    trees$max.samp %<>% signif(digits = 2)
    
    trees
  }
  
  samp.trees
  
}
random.crown <- function(tree.clips, tree.df, cores) {
  # Accepts a list of file names (each one corresponding to a laz file
  # centered on a randomly selected tree), and returns a data frame containing crown variables
  # for each randomly selected height associated with the focal tree. Can be parallelized
  # by supplying the number of cores over which computations can be spread.
  
  require(foreach)
  require(doMC)
  registerDoMC(cores = cores)
  
  set.seed(666)
  stats <- foreach(i = 1:nrow(tree.df), .combine = rbind, .errorhandling = 'remove') %dopar% {
    
    xy <- data.frame(x = tree.df[i, 'coords.x1'], y = tree.df[i, 'coords.x2'])
    
    las <- readLAS(tree.clips[[i]])
  
    # Height-sampling is set to select from a height of 10 m, to the maximum sampling height
    # identified for each tree. Change "10" to a different value if you'd like to set it to something
    # different. For now, it is hard-coded in the function.
    ht.samp <- sample(10:tree.df[i, 'max.samp'], size = 10)
    hts <- data.frame(tree = i,
                      height = ht.samp,
                      exp.top = tree.df[i, 'Z'])
    hts$tree %<>% as.character()
    hts$height %<>% as.numeric() %>% sort()
    hts$closure <- 0
    
    ground <- filter_poi(las, Classification == 2)
    # grid_terrain arguments have been found to work satisfactorily, but may need to be
    # adjusted to work with other datasets.
    dtm <- grid_terrain(ground, res = 0.1, algorithm = knnidw(k = 10))
    elev <- raster::extract(dtm, xy, 'bilinear')
    
    slope.rast <- terrain(dtm, opt = 'slope')
    slope <- slope.rast %>% raster::as.data.frame(na.rm = T) %>% .$slope %>% mean(na.rm = T)
    
    # Although we have the absolute tree top xyz position from the initial tree detection approach,
    # we also use the same "top-finding" heuristic to find the tree top.
    # I think this approach is needless, though--why not just use the exact coordinates for both
    # random and study tree tops?
    veg <- clip_circle(las, xcenter = tree.df[i, 1], ycenter = tree.df[i, 2], radius = 0.1) %>%
      filter_poi(Classification == 1)
    max <- veg@data$Z %>% max()
    top <- max - elev
    
    norm <- normalize_height(las, knnidw()) %>% filter_poi(Classification == 1)
    zmed <- norm@data$Z %>% median()
    
    las %<>% filter_poi(Classification == 1 & Z > (10 + elev))
    
    rho.in <- lidR::density(las)
    adj.rho.in <- rho.in / zmed
    
    n.points <- npoints(las)
    
    for(j in 1:nrow(hts)){
      hts[j, 'closure'] <- filter_poi(las, Z > (elev + hts[j, 2])) %>%
        npoints()
    }
    
    hts$elevation <- elev
    hts$slope <- slope
    hts$top <- top
    hts$depth <- hts$top - hts$height
    hts$top.diff <- hts$top - hts$exp.top
    hts$zmed <- zmed
    hts$rho.in <- rho.in
    hts$adj.rho.in <- adj.rho.in
    hts$n.points <- n.points
    
    hts
  }
  
  stats
  
}
study.crown <- function(tree.clips, trees.df, hts.df, cores) {
  # This function accepts a list of file names (each one corresponding to a laz file
  # centered on a climbed tree), a data frame that associates a tree with an xy position,
  # and another data frame that describes all of the sample height within a tree. It returns a data frame containing crown variables
  # for each sample height associated with the focal tree. Operations can be parallelized to run over
  # a specified number of cores, but this is excessive for eight trees.
  
  require(foreach)
  require(doMC)
  registerDoMC(cores = cores)
  
  set.seed(666)
  stats <- foreach(i = 1:nrow(trees.df), .combine = rbind) %dopar% {
    
    xy <- data.frame(x = trees.df[i, 2], y = trees.df[i, 3])
    
    las <- readLAS(tree.clips[[i]])
    
    tree <- trees.df[i, 1]
    hts <- hts.df[hts.df$tree == tree, ]
    
    hts$tree %<>% as.factor()
    hts$height %<>% as.numeric() %>% sort()
    hts$closure <- 0
    
    ground <- filter_poi(las, Classification == 2)
    dtm <- grid_terrain(ground, res = 0.5, algorithm = knnidw(k = 10))
    elev <- raster::extract(dtm, xy, 'bilinear')
    
    slope.rast <- terrain(dtm, opt = 'slope')
    slope <- slope.rast %>% raster::as.data.frame(na.rm = T) %>% .$slope %>% mean(na.rm = T)
    
    # Although we have the absolute tree top xyz position from the initial tree detection approach,
    # we also use the same "top-finding" heuristic to find the tree top.
    # I think this approach is needless, though--why not just use the exact coordinates for both
    # random and study tree tops?
    veg <- clip_circle(las, xcenter = trees.df[i, 2], ycenter = trees.df[i, 3], radius = 0.5) %>%
      filter_poi(Classification == 1)
    max <- veg@data$Z %>% max()
    top <- max - elev
    
    norm <- normalize_height(las, knnidw()) %>% filter_poi(Classification == 1)
    zmed <- norm@data$Z %>% median()
    
    rho.in <- lidR::density(las)
    adj.rho.in <- rho.in / zmed
    
    n.points <- npoints(las)
    
    for(j in 1:nrow(hts)){
      
      ht <- hts[j, 4]
      
      hts[j, 'closure'] <- filter_poi(las, Z > (elev + ht)) %>%
        npoints() %>% divide_by(n.points) %>% multiply_by(zmed)
    }
    
    hts$x <- xy$x
    hts$y <- xy$y
    hts$elevation <- elev
    hts$z <- elev + hts$height
    hts$slope <- slope
    hts$top <- top
    hts$depth <- top - hts$height
    hts$zmed <- zmed
    hts$rho.in <- rho.in
    hts$adj.rho.in <- adj.rho.in
    hts$n.points <- n.points
    
    las@data$Z <- las@data$Z - elev
    file.name <- paste0(tree, '.adj.laz')
    writeLAS(las, here::here(out.path, 'clips.adj', file.name))
    
    hts
  }
  
  stats
  
}
loess.index <- function(all.crowns.df, form, span = 0.5, degree = 1, df.out = T){
  # This function accepts a data frame containing the height, depth, and closure
  # values for either virtually sampled or actually sampled trees, fits a LOESS
  # curve for a supplied span and degree, and optionally outputs a data frame
  # or the raw fit as output.
  
  form %<>% as.formula()
  
  fit <- loess(form, all.crowns.df,
               span = span, degree = degree
  )
  if(df.out == T){
    fit.df <- data.frame(fitted = fit$fitted, residuals = fit$residuals)
    cbind(all.crowns.df, fit.df) 
  }
  else if(df.out == F){
    fit
  }
}
kmean.it <- function(all.crowns.df, var, k = 2, iterate = 999, starts = 999){
  stats <- data.frame(var = all.crowns.df[, var])
  
  set.seed(666)
  stats %<>% kmeans(centers = k, iter.max = iterate, nstart = starts)
  
  all.crowns.df$group <- stats$cluster
  
  all.crowns.df
  
}

# Phyloseq processing ####
# See use in: "compile.R", clean.R" and "composition.R"
remove.contam <- function(phy, controls, prop) {
  # Remove otus with >1% of experimental reads found in the pcr and extraction controls
  # This function takes a phyloseq object, a list of OTUs found in controls, and an acceptable proportion of reads
  # for a given OTU (among all reads in the dataset) that can be found in the controls. OTUs exceeding this 
  # proportion are removed from the input phyloseq, and this change is reflected in the output phyloseq object.
  
  # Get the current named list of otu read totals
  otu.tab <- phy@otu_table %>% data.frame()
  otu.names <- colnames(otu.tab)
  
  # Identify otus found in the negative controls
  contam.sums <- otu.tab[, colnames(otu.tab) %in% names(controls)] %>% colSums()
  
  # Calculate a ratio of the number of reads (in negative controls):(in non-control samples) for each otu
  controls <- controls[names(controls) %in% names(contam.sums)]
  
  contam.prop <- controls / contam.sums
  
  # Get the names of the otus with ratio > prop and remove them from the input phyloseq object
  contam.names <- contam.prop[contam.prop > prop] %>% names()
  
  keep <- otu.names[!(otu.names %in% contam.names)]
  
  phy %<>% prune_taxa(keep, .)
  
  keep <- sample_sums(phy) > 0
  prune_samples(keep, phy)
  
}
perfect.filt <- function(phy, k, file.prefix) {
  # Taking a phyloseq oject as input, this script removes OTUs which make non-significant contributions (alpha = 0.1) to the total
  # covariance of the dataset, producing a filtered phyloseq object and a summary rds + figure as outputs.
  
  # For reproducible filtering
  set.seed(666)
  
  otu.tab <- phy@otu_table %>% data.frame()
  
  # Even though we will ultimately use permutation filtering, the authors mention that performance is improved with Order = 'pvals'.
  # This requires that we go ahead and obtain simulataneous PERFect output first. 
  sim <- PERFect_sim(X = otu.tab)
  # Of note here: by default, an alpha of 0.1 is used. Taxa with p-values greater than 0.1 are filtered out.
  # For this reason, filtering is somewhat conservative, which would ideally allow us to retain real and somewhat rare taxa.
  perm <- PERFect_perm(X = otu.tab, Order = 'pvals', pvals_sim = sim, algorithm = 'full', k = k)
  tab.out <- perm$filtX
  
  # Save PERFect output
  here(out.path, paste0(file.prefix, '.', 'perfect.perm.rds')) %>% saveRDS(perm, .)
  
  # Make p-value plots and save to output
  perm %<>% pvals_Plots(otu.tab)
  perm <- perm$plot + ggtitle(paste0('Permutation filtering', '-', file.prefix)) + scale_color_colorblind()
  here(out.path, paste0(file.prefix, '.', 'perfect.perm.pvals.pdf')) %>% ggsave(., perm, dpi = 300)
  
  # Update the phyloseq object
  keep <- tab.out %>% colnames()
  phy %<>% prune_taxa(keep, .)
  otu_table(phy) <- otu_table(tab.out, taxa_are_rows = F)
  keep <- sample_sums(phy) > 0
  prune_samples(keep, phy)
  
}
parse.tax <- function(phy) {
  # Taking a phyloseq object as input, this script coverts the UNITE taxonomy to something more readable.
  
  phy@tax_table@.Data %<>% parse_taxonomy_greengenes()
  
  tax.tab <- phy@tax_table %>% data.frame()
  tax.tab$Genus_species <- ifelse(!is.na(tax.tab$Species), paste0(tax.tab$Genus, ' ', tax.tab$Species),
                                  ifelse(!is.na(tax.tab$Genus), paste0(tax.tab$Genus, ' sp.'),
                                         ifelse(!is.na(tax.tab$Family), paste0(tax.tab$Family, ' sp.'),
                                                ifelse(!is.na(tax.tab$Order), paste0(tax.tab$Order, ' sp.'),
                                                       ifelse(!is.na(tax.tab$Class), paste0(tax.tab$Class, ' sp.'),
                                                              ifelse(!is.na(tax.tab$Phylum), paste0(tax.tab$Phylum, ' sp.'),                                                                      ifelse(!is.na(tax.tab$Kingdom), paste0(tax.tab$Kingdom, ' sp.'), NA)))))))
  phy@tax_table <- tax.tab %>% as.matrix() %>% tax_table()
  
  phy
  
}
filt.n.prune <- function(phy.in, prop = 0, depth = 0) {
  # This function accepts a phyloseq object as input, removes OTUs that are not
  # present in a specified proportion of samples (prop), and removes samples
  # that fail to exceed a specified depth afterwards. A filtered phyloseq object is returned.
  
  # Care should be taken to determine whether filt.n.prune() is being applied to relative
  # abundance or count data! I've generally only set prop and depth to zero where it's called,
  # so it's not always apparent from composition.R how this function works!
  filt.phy <- filter_taxa(phy.in, function(x) {sum(x > 0) > (prop*length(x))}, prune = T)
  keep <- sample_sums(filt.phy) > depth
  filt.phy %>% prune_samples(keep, .)
  
}
joiner <- function(phy.in){
  # This function take a phyloseq object as input, makes some aesthetic changes
  # related to variable coding, variable inclusion, and factor level character
  # strings. It has a very specific use case and will ideally be removed in future
  # applications, as its only purpose is to correct earlier data formatting choices.
  
  sample_data(phy.in)$height <- sample_data(phy.in)$ht
  meta.tab <- phy.in@sam_data %>% data.frame() %>%
    dplyr::select(sampleID, age, height, tree, tree.ht)
  meta.tab$tree %<>% str_replace('DT_NEIGHBOR', 'DT NEIGHBOR')
  meta.tab %<>% left_join(k.crowns, by = c('tree', 'height'))
  rownames(meta.tab) <- meta.tab$sampleID
  sample_data(phy.in) <- sample_data(meta.tab)
  phy.in
}
phy.clean <- function(phy.in, depth = 0){
  # This function accepts a phyloseq object as input and returns a list of phyloseq objects
  # that have received different levels of processing.
  
  counts <- age.split(phy.in, depth = depth)
  relative <- lapply(counts, transform_sample_counts, fun = rel)
  logged <- lapply(relative, gen.log)
  # # Apply group.split only works like this because none of the filtering operates at the level of reads in a given sample.
  # # In theory, age.split, group.split, and tree.split can happen anywhere in the pipeline!
  # grouped <- lapply(logged, group.split)
  # tree.log <- tree.split(counts$all) %>% lapply(transform_sample_counts, fun = rel) %>% lapply(gen.log)
  # list(counts = counts, relative = relative, logged = logged,
  #      grouped = grouped, tree.log = tree.log)
  list(counts = counts, relative = relative, logged = logged)
}
age.split <- function(phy.in, prop = 0, depth = 0){
  # This function takes a phyloseq object as input, removes sample heights that failed to occur in the first age class,
  # separates the input according needle age class, and returns a list of filtered phyloseq objects associated
  # with the full dataset and each age class subset individually.
  
  all <- phy.in %>% subset_samples(tree.ht %in% a1.tree.ht) %>% filt.n.prune(prop = prop, depth = depth)
  
  a1 <- phy.in %>% subset_samples(age == 'A1') %>% subset_samples(tree.ht %in% a1.tree.ht) %>%
    filt.n.prune(prop = prop, depth = depth)
  a2 <- phy.in %>% subset_samples(age == 'A2') %>% subset_samples(tree.ht %in% a1.tree.ht) %>%
    filt.n.prune(prop = prop, depth = depth)
  a3 <- phy.in %>% subset_samples(age == 'A3') %>% subset_samples(tree.ht %in% a1.tree.ht) %>%
    filt.n.prune(prop = prop, depth = depth)
  a4 <- phy.in %>% subset_samples(age == 'A4') %>% subset_samples(tree.ht %in% a1.tree.ht) %>%
    filt.n.prune(prop = prop, depth = depth)
  list(all = all, a1 = a1, a2 = a2, a3 = a3, a4 = a4)
}
tree.split <- function(phy.in, prop = 0, depth = 0){
  # This function receives a phyloseq object as input, removes sample heights that failed to occur in the first age class,
  # separates the input according tree, and returns a list of filtered phyloseq objects associated
  # with the full dataset and each tree subset individually.
  
  all <- phy.in %>% subset_samples(tree.ht %in% a1.tree.ht) %>% filt.n.prune(prop = prop, depth = depth)
  phy.in %<>% subset_samples(age != 'A1') %>% subset_samples(tree.ht %in% a1.tree.ht)
  BIRD348 <- phy.in %>% subset_samples(tree == 'BIRD348') %>%
    filt.n.prune(prop = prop, depth = depth)
  CC <- phy.in %>% subset_samples(tree == 'CC') %>%
    filt.n.prune(prop = prop, depth = depth)
  DISCOVERY <- phy.in %>% subset_samples(tree == 'DISCOVERY') %>%
    filt.n.prune(prop = prop, depth = depth)
  DT.NEIGHBOR <- phy.in %>% subset_samples(tree == 'DT NEIGHBOR') %>%
    filt.n.prune(prop = prop, depth = depth)
  LATREEN <- phy.in %>% subset_samples(tree == 'LATREEN') %>%
    filt.n.prune(prop = prop, depth = depth)
  PC02 <- phy.in %>% subset_samples(tree == 'PC02') %>%
    filt.n.prune(prop = prop, depth = depth)
  PC12 <- phy.in %>% subset_samples(tree == 'PC12') %>%
    filt.n.prune(prop = prop, depth = depth)
  PC17 <- phy.in %>% subset_samples(tree == 'PC17') %>%
    filt.n.prune(prop = prop, depth = depth)
  list(all = all,
       BIRD348 = BIRD348, CC = CC, DISCOVERY = DISCOVERY, DT.NEIGHBOR = DT.NEIGHBOR,
       LATREEN = LATREEN, PC02 = PC02, PC12 = PC12, PC17 = PC17)
}
group.split <- function(phy.in, prop = 0, depth = 0){
  # This function receives a phyloseq object as input, separates the input according exposure group,
  # and returns a list of filtered phyloseq objects associated
  # with the full dataset and each exposure group subset individually.
  
  open <- phy.in %>% subset_samples(group == 'Open') %>% filt.n.prune(prop = prop, depth = depth)
  closed <- phy.in %>% subset_samples(group == 'Closed') %>% filt.n.prune(prop = prop, depth = depth)
  list(Open = open, Closed = closed)
}
# Transformations
rel <- function(x) {
  # Transform count data to relative abundance
  
  x/(sum(x))
  
}
pa <- function(x) {
  
  # Transform quantitative data to presence-absence data.
  
  ifelse(x > 1, 1, 0)
  
}
gen.log <- function(phy.in){
  # Adapted form McCune and Grace's "Analysis of Ecological Communities",
  # this function accepts a phyloseq object as input, performs a generalized
  # log transform (base e), and returns a transformed phyloseq object as output.
  otu.tab <- phy.in@otu_table %>% data.frame()
  otu.tab[otu.tab == 0] <- NA
  c <- min(otu.tab, na.rm = T) %>% log() %>% as.integer()
  d <- exp(c)
  transform_sample_counts(phy.in, fun = function(x){log(x + d) - c})
}

# Taxa ####
# See use in: "clean.R" and "taxa.R"
rank.n.stack <- function(phy.in, var, prop){
  
  age <- phy.in@sam_data %>% data.frame() %>% .$age %>% unique()
  if(length(age) > 1){
    age <- 'All'
  }
  
  tax.table <- phy.in@tax_table %>% data.frame() %>%
    rownames_to_column(var = 'OTU')
  
  counts.merged <- phy.in %>%
    merge_samples(var) %>%
    psmelt() %>% select(OTU, Sample, Abundance)
  
  rel.merged <- phy.in %>%
    merge_samples(var) %>%
    transform_sample_counts(function(x) {x/sum(x)}) %>%
    psmelt() %>% select(OTU, Sample, Abundance)
  
  name <- rel.merged$Sample %>% unique() 
  
  out <- NULL %>% as.list()
  for(i in 1:length(name)){
    
    var.level <- filter(rel.merged, Sample == name[[i]])
    var.level <- var.level[order(var.level$Abundance, decreasing = T), ]
    var.level$total <- var.level$Abundance
    
    for(j in 1:length(var.level$OTU)) {
      if(j == 1) {
        var.level[1, 4] <- var.level[1, 3]
      }
      else {
        var.level[j, 4] <- var.level[j, 3] + var.level[(j-1), 4]
      }
    }
    
    var.level %<>% filter(total <= prop)
    
    var.level$Abundance <- NULL
    var.level$total <- NULL
    
    var.level %<>% left_join(counts.merged, by = c('OTU', 'Sample'))
    reads <- var.level %>% group_by(Sample) %>% summarize(total = sum(Abundance))
    var.level %<>% left_join(reads, by = 'Sample')
    var.level$Abundance <- var.level$Abundance / var.level$total
    
    out[[i]] <- var.level
    
  }
  
  out %<>% bind_rows() %>% left_join(tax.table, by = 'OTU')
  out$Sample %<>% str_replace('DT_NEIGHBOR', 'DT NEIGHBOR')
  out$Genus %<>% str_replace('Rasutoria', 'Zasmidium')
  out$Order[out$Genus == 'Micraspis'] <- 'Micraspidales'
  
  out$combo <- glue::glue('<i>{out$Genus}</i> | {out$Order}')
  out[out == '<i>NA</i> | NA'] <- 'NA'
  out$combo %<>% str_replace('<i>NA</i>', 'NA')
  
  out <- out[order(out$Order), ]
  out$age <- age
  order <- out$combo %>% unique()
  out$combo <- factor(out$combo, levels = c(order))
  
  out
  
}
sad.occ.abund <- function(phy.in, var){
  # This function accepts a phyloseq object as input, merges samples according to
  # a user supplied "var" argument (must be categorical), and returns a two-element
  # list of plot-ready data frames. The first data frame describes the
  # species abundance distribution ("sad") across this var, and the second describes
  # the occupancy-abundance relationship using all OTUs ("otus").
  
  merged <- phy.in %>%
    merge_samples(var)
  
  otus.pa <- merged %>%
    transform_sample_counts(fun = pa) %>%
    otu_table() %>%
    data.frame() %>% colSums()
  
  otus.pa
  
  sad <- otus.pa %>% as.character() %>% data.frame(samples = .) %>%
    group_by(samples) %>% summarize(otus = n())
  sad$samples %<>% as.numeric()
  sad$otus <- (sad$otus / sum(sad$otus)) %>% round(digits = 2)
  
  otus.ra <- phy.in %>%
    transform_sample_counts(fun = rel) %>%
    otu_table() %>%
    data.frame() %>% colMeans()# %>% log()
  
  otus <- data.frame(pa = otus.pa, ra = otus.ra)
  
  list(sad = sad, otus = otus)
  
}
rare.curve <- function(phy, step = 10){
  
  set.seed(666)
  
  otu.tab <- phy@otu_table %>% data.frame()
  rarecurve(otu.tab, step=step, label = F, xlim = c(0, 5000), alpha = 0.5)
  
  rare.tab <- lapply(1:length(rare), FUN = function(x){
    cbind(rep(x, length(rare[[x]])), attr(rare[[x]], 'Subsample'),
          rare[[x]]) %>%
      as.data.frame(row.names = NULL)
  }) %>% do.call('rbind', .)
  
  rare.tab %<>% dplyr::select(sample = V1,
                              depth = V2,
                              ASVs = V3)
  
  rare.tab$sample %<>% as.factor()
  
  ggplot(rare.tab, aes(x = depth, y = ASVs, color = sample)) +
    geom_line(alpha = 0.5, color = 'black') +
    theme(legend.position = 'none')
  
} # see if still relevant/needed

# Diversity ####
# As used in "diversity.R"
inext.div <- function(in.phy, level = 1000, cores){
  # This function accepts a phyloseq object as input, estimates the diversity of each
  # sample for Hill numbers 0-2 at user-defined sampling depth (level), and returns
  # a data frame associated diversity calculations with existing sample metadata. This
  # function can be parallelized to operate across a specified number of cores.
  
  require(foreach)
  require(doMC)
  registerDoMC(cores = cores)
  
  otu.tab <- in.phy@otu_table %>% t %>% data.frame()
  meta.tab <- in.phy@sam_data %>% data.frame() %>% rownames_to_column(var = 'sample')
  meta.tab$sample %<>% str_replace('-', '.')
  
  div <- foreach(i = 1:ncol(otu.tab), .combine = bind_rows) %dopar% {
    set.seed(666)
    estimateD(otu.tab[, i], level = level) %>%
      tibble() %>%
      mutate(sample = colnames(otu.tab)[i])
  }
  
  div %<>%
    transmute(sample = sample,
              metric = case_when(order == 0 ~ 'q0',
                                 order == 1 ~ 'q1',
                                 order == 2 ~ 'q2'),
              estimate = qD,
              lcl = qD.LCL,
              ucl = qD.UCL)
  
  # It's not apparent to me how we can prevent q2 from being calculated in the first place...
  div %>% left_join(meta.tab, by = 'sample') %>% filter(metric != 'q2')
  
}
anova.crown <- function(x, form, adjust = F){
  # This function accepts a grouped data frame and an ANOVA model formula as input,
  # performs an ANOVA testing the supplied model for estimated richness and Shannon index,
  # cleans up the test output, and returns a results data frame as the final output.
  
  metric <- x$metric %>% unique()
  if(metric == 'q0'){
    title <- 'Estimated richness'
  }
  else if(metric == 'q1'){
    title <- 'Estimated Shannon index'
    x$estimate %<>% log2()
  }
  
  form %<>% as.formula()
  model <- aov(form, data = x)
  test <- car::Anova(model, white.adjust = adjust)
  if(adjust == T){
    
    test %<>% data.frame() %>% rownames_to_column(var = 'Term') %>%
      select(Term, Df, `F`, P = Pr..F.)
    
  } else {
    
    test %<>% data.frame() %>% rownames_to_column(var = 'Term') %>%
      select(Term, Df, `F` = F.value, P = Pr..F.)
    
  }

  test[,3:4] <- round(test[, 3:4], digits = 3)

  test$P %<>% as.character()
  test[test < 0.001] <- '< 0.001'

  total <- data.frame(Term = 'Total', Df = sum(test$Df))
  test %<>% bind_rows(total)
  
  test$metric <- title
  
  test
}
tukey <- function(x, form) {
  # This function accepts a grouped data frame and a single-term model formula as input,
  # perform Tukey Honestly Significant Differences tests for all pairwise comparisons
  # between all levels of the supplied model term, assigns the level a significance letter,
  # and outputs a data frame that can be incorporated into plots.
  
  require(multcompView)
  
  metric <- x$metric %>% unique()
  if(metric == 'q1'){
    x$estimate %<>% log2()
  }
  
  var <- form
  form <- paste0('estimate ~ ', form) %>% as.formula()
  
  model <- aov(form, data = x)
  levs <- model$xlevels[[var]]
  hsd <- TukeyHSD(model)
  
  p <- hsd[[var]][, 4]
  labs <- data.frame(multcompLetters(p)['Letters'], stringsAsFactors = F)
  labs %<>% rownames_to_column(var)
  labs %<>% .[match(levs, .[[var]]) ,]
  labs[[var]] %<>% factor(levels = levs)
  
  letters <- labs$Letters %>% as.character()%>%
    strsplit('') %>% unlist %>% unique
  
  for(i in 1:nrow(labs)){
    foo <- labs$Letters[i] %>% as.character() %>% strsplit('') %>% unlist
    labs$Letters[i] <- letters[which(letters %in% foo)] %>% paste0(collapse = '')
  }
  
  labs %>% tibble(metric = unique(x$metric))
}
mean.boot <- function(x, y, var, n.perm = 999){
  # This function accepts a grouped data frame, a categorical variable of interest (var),
  # and the number of bootstrap permutations (n.perm) as inputs, calculates bootstrapped mean
  # alpha diversity values, and returns tibble columns of these values as output.
  
  metric <- y$metric %>% unique()
  if(metric == 'q1'){
    x$estimate %<>% log2()
  }
  
  boots <- replicate(n.perm, mean(sample(x$estimate, replace = T)))
  
  tibble(metric = y$metric,
         factor = y[, var, drop = T],
         estimate = boots)
}

# Composition ####
# As used in "composition.R"
# Ordination
nmds.mc.par <- function(phy.in, k = 2, distance = 'bray', trymax = 50, autotransform = FALSE, no.share = F,
                        trace = 1, zerodist = 'ignore', display = 'sites', test = T, perm = 999, col.hist = 'blue',
                        col.line = 'red', lty = 2, las = 1, lab = c(5, 5, 4), n.cores, ...){
  
  # This function accepts a phyloseq object formatting preferences, and several vegan::metaMDS() arguments as inputs,
  # performs non-metric multidimensional scaling with metaMDS, performs a Monte-Carlo permutation test with a user-supplied
  # number of permutations (similar to PC-ORD), and returns MC summary plots and list of ordination and testing
  # result objects as output. The Monte Carlo takes a while to compute, and can be toggled on or off (test = T or F).
  # Operation of the function can also be parallelized over a specified number of cores.

  require(vegan)
  require(MASS)
  require(foreach)
  require(doMC)
  registerDoMC(cores = n.cores)
  
  x <- phy.in@otu_table %>% data.frame()
  meta.tab <- phy.in@sam_data %>% data.frame() %>% rownames_to_column(var = 'sample')
  set.seed(666)
  z <- metaMDS(comm = x, k = k, distance = distance, trymax = trymax,
             autotransform = autotransform, noshare = no.share, trace = trace, zerodist = zerodist, ...) #nmds analysis
  
  # Get stress
  z.stress <- z$stress
  
  if(test == T){
    
    # Permutation test that NMDS stress does not differ from randomly simulated stress
    # Obtain random final stresses from n = perm sets of column-shuffled data
    y.stress <- foreach(i = 1:perm, .combine = c) %dopar% {
      apply(x, 2, sample) %>%
        metaMDS(k = k, distance = distance, trymax = trymax,
                autotransform = autotransform, noshare = no.share, trace = 0, zerodist = zerodist, ...) %>%
        .$stress
    }
    
    # Compute number of random runs with stress <= observed
    n <- sum(y.stress <= z.stress)
    # Compute p-value
    p.value <- (1 + n) / (1 + perm)
    
    # Make histogram of simulated stresses
    xmin <- min(z.stress, min(y.stress))
    xmax <- max(z.stress, max(y.stress))
    hist(y.stress, col = col.hist, las = las, lab = lab,
         xaxs = 'i', yaxs = 'i', xlim = c(xmin,xmax), xlab = 'Stress',
         main = paste('Random Permutation Distribution of Stress for',
                      k, 'Dimensions', sep = ' '), ...)
    abline(v = z.stress, col = col.line, lty = lty, lwd = 2, ...)
    
    # Print out results of the permutation test 
    cat('Randomization Test of Stress:\n')
    cat('Permutation stress values:\n')
    print(y.stress)
    rbind('Observed stress' = z.stress, 'P-value' = p.value) %>% print()
    
    p.value %<>% signif(digits = 3)
    
  }
  
  else if(test == F){
    
    p.value <- 'NA'
    
  }
  
  z.stress %<>% signif(digits = 3)
  
  # Output NMDS scores
  scores <- scores(z, display = display) %>%
    data.frame() %>% rownames_to_column(var = 'sample') %>%
    left_join(meta.tab, by = 'sample')
  
  list(ord = z, scores = scores, stress = z.stress, p.value = p.value)
}
nmds.plot <- function(nmds.in, note = T, ...){
  # This function accepts an NMDS ordination object (produced in vegan) as input, returning a
  # skeletal ggplot2 object (amenable to customization in the scripting environment) as output. The results of the
  # Monte-Carlo permutation test performed by nmds.mc.par() can be included or omitted (note = T or F)
  
  scores <- nmds.in$scores
  stress <- nmds.in$stress
  nmds.pv <- nmds.in$p.value
  
  axes <- list(...)
  
  if('NMDS1' %in% axes & 'NMDS2' %in% axes){
    
    plot <- ggplot(scores, aes(x = NMDS1, y = NMDS2)) +
      xlab('\nNMDS1') +
      ylab('NMDS2\n')
    
  }
  else if('NMDS1' %in% axes & 'NMDS3' %in% axes){
    scores$NMDS2 <- NULL
    
    plot <- ggplot(scores, aes(x = NMDS1, y = NMDS3)) +
      xlab('\nNMDS1') +
      ylab('NMDS3\n')
    
  }
  
  else if('NMDS2' %in% axes & 'NMDS3' %in% axes){
    scores$NMDS1 <- NULL
    
    plot <- ggplot(scores, aes(x = NMDS2, y = NMDS3)) +
      xlab('\nNMDS2') +
      ylab('NMDS3\n')
    
  }
  
  nmds.lab <- paste0('NMDS: Stress = ', stress, ', P = ', nmds.pv)
  
  if(note == T){
    plot <- plot + annotate('text', x = min(scores[, 2]), y = min(scores[, 3]),
                            label = nmds.lab, hjust = 0, size = 2.5)
  }
  
  plot +
    theme_cowplot() +
    theme(axis.title.x = element_text(size = 7, face = 'bold'),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 7, face = 'bold'),
          axis.text.y = element_text(size = 7),
          legend.title = element_text(size = 7, face = 'bold'),
          legend.text = element_text(size = 7),
          plot.background = element_rect(fill = 'white', color = 'white')
          )
}
dbrda.it <- function(log.in, rel.in = NULL, form, dist = 'bray', vectors = F, min.r = NULL, max.p = NULL){
  # This function accepts a phyloseq object, an optional relative abundance-transformed phyloseq object for taxonomic vector
  # creation, a distance/dissimilarity calculation method (dist), and minimum acceptable r2 and maximum p-values for vectors.
  # Vectors can be toggled on or off (vectors = T or F). A skeletal ggplot2 object (amenable to customization in the scripting
  # environment) is returned as output.
  
  otu.tab <- log.in@otu_table %>% data.frame()
  # meta.tab <- log.in@sam_data %>% data.frame() %>% dplyr::select(sampleID, tree, age, height, depth, closure, residuals, group)
  # meta.tab <- log.in@sam_data %>% data.frame() %>% dplyr::select(sampleID, tree, age, height, depth, closure, group)
  meta.tab <- log.in@sam_data %>% data.frame() %>% dplyr::select(sampleID, tree, age, height, closure, group)
  # meta.tab %<>% cbind(otu.tab)
  form %<>% as.formula()
  
  obj <- dbrda(formula = form, distance = dist, data = meta.tab)
  stats <- summary(obj)
  
  prop <- stats$cont[[1]] %>% data.frame()
  db.1.prop <- prop$dbRDA1[2] %>% round(digits = 3) * 100
  db.2.prop <- prop$dbRDA2[2] %>% round(digits = 3) * 100
  # db.2.prop <- prop$MDS1[2] %>% round(digits = 3) * 100
  
  db.1.title <- paste0('\ndbRDA1 [', db.1.prop, '%]')
  db.2.title <- paste0('dbRDA2 [', db.2.prop, '%]\n')
  # db.2.title <- paste0('MDS1 [', db.2.prop, '%]\n')
  
  sites <- scores(obj, display = 'sites') %>% data.frame() %>%
    rownames_to_column(var = 'sampleID') %>% left_join(meta.tab, by = 'sampleID') %>%
    column_to_rownames(var = 'sampleID')
  
  if(vectors == F){
    plot <- ggplot(data = sites, mapping = aes(x = dbRDA1, y = dbRDA2)) +
      # plot <- ggplot(data = sites, mapping = aes(x = dbRDA1, y = MDS1)) +
      xlab(db.1.title) +
      ylab(db.2.title) +
      theme_cowplot() +
      theme(axis.title.x = element_text(size = 7, face = 'bold'),
            axis.text.x = element_text(size = 7),
            axis.title.y = element_text(size = 7, face = 'bold'),
            axis.text.y = element_text(size = 7),
            legend.title = element_text(size = 7, face = 'bold'),
            legend.text = element_text(size = 7),
            legend.text.align = 0.5)
  }
  
  else if(vectors == T){
    
    vect <- stats$biplot %>% data.frame() %>% rownames_to_column(var = 'var')
    set.seed(666)
    term.test <- anova(obj, by = 'margin', strata = meta.tab$tree) %>% data.frame()
    vect$p <- term.test[1:2, 4]
    # vect$p <- term.test[1, 4]
    vect$r <- 0.1 # Numeric placeholder
    vect %<>% dplyr::select(-MDS4, -MDS1, -MDS2, -MDS3)
    # vect %<>% dplyr::select(-MDS4, -MDS5, -MDS2, -MDS3)
    
    vect.otu.tab <- rel.in@otu_table %>% data.frame() %>% rownames_to_column(var = 'sampleID') %>% left_join(meta.tab, by = 'sampleID')
    # vect.otu.tab %<>% dplyr::select(height, closure, otu.1, otu.2, otu.3, otu.6)
    vect.otu.tab %<>% dplyr::select(otu.1, otu.2, otu.3, otu.6)
    
    set.seed(666)
    fit <- envfit(obj, vect.otu.tab, strata = meta.tab$tree)
    # adj.p.val <- p.adjust(fit$vectors$pvals, method = 'fdr')
    # fit$vectors$pvals <- adj.p.val
    
    # fit.filt <- data.frame(r = fit$vectors$r, p = fit$vectors$pvals) %>%
    #   rownames_to_column(var = 'var') %>%
    #   filter(r > min.r, p < max.p) %>% dplyr::select(-r) %>% left_join(tax.tab, by = 'var')
    
    fit.filt <- data.frame(p = fit$vectors$pvals, r = fit$vectors$r) %>%
      rownames_to_column(var = 'var')
    
    fit.scores <- as.data.frame(scores(fit, display = 'vectors')) %>%
      rownames_to_column(var = 'var') %>% right_join(fit.filt, by = 'var')
    
    # sig <- fit.scores %>% filter(p < max.p, r > min.r) %>% dplyr::select(-r) %>% 
    #   rbind(vect) %>% filter(p < max.p)
    # fit.scores$r <- NULL
    
    vect %<>% rbind(fit.scores)
    # vect <- fit.scores
    sig <- vect %>% filter(p < max.p, r > min.r)
    
    plot <- ggplot(data = sites, mapping = aes(x = dbRDA1, y = dbRDA2)) +
    # plot <- ggplot(data = sites, mapping = aes(x = dbRDA1, y = MDS1)) +
    geom_segment(data = vect,
                 aes(x = 0,
                     xend = dbRDA1 * 2,
                     y = 0,
                     yend = dbRDA2 * 2,
                     # yend = MDS1 * 2,
                     color = var),
                 alpha = 0.7,
                 arrow = arrow(length = unit(2, 'mm')),
                 size = 0.75) +
    geom_point(aes(x = dbRDA1 * 2.2,
                   y = dbRDA2 * 2.2
                   # y = MDS1 * 2.2
                   ),
               sig,
               shape = 8,
               size = 1) +
    xlab(db.1.title) +
    ylab(db.2.title) +
    theme_cowplot() +
    theme(axis.title.x = element_text(size = 7, face = 'bold'),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 7, face = 'bold'),
          axis.text.y = element_text(size = 7),
          legend.title = element_text(size = 7, face = 'bold'),
          legend.text = element_text(size = 7),
          legend.text.align = 0)
  }
  
  return(plot)
  
}

# Analyses
permdisp <- function(phy.in, test = 'group', method = 'bray', type = 'median', n.perm = 999, bias = F, pair = F, block = 'none') {
  # This function accepts phyloseq objects as inputs, performs PERMDISP,
  # and returns formatted data frames for visualization.
  # Relevant phyloseq::distance, vegan::betadisper, and vegan::permutest options can be set as well.
  
  set.seed(666)
  
  d <- distance(phy.in, method = method)
  meta.tab <- phy.in@sam_data %>% data.frame()
  
  if(test == 'group') {
    group <- phy.in@sam_data$group
  }
  else if(test == 'tree') {
    group <- phy.in@sam_data$tree
  }
  else if(test == 'age') {
    group <- phy.in@sam_data$age
  }
  
  perm <- how(nperm = n.perm)
  if(block == 'tree'){
    setBlocks(perm) <- with(meta.tab, tree) 
  }
  
  else if(block == 'group'){
    setBlocks(perm) <- with(meta.tab, group) 
  }
  
  else if(block == 'age'){
    setBlocks(perm) <- with(meta.tab, age) 
  }
  
  bd <- betadisper(d = d, group = group, type = type, bias.adjust = bias)
  df <- permutest(bd, pairwise = pair, permutations = perm)$tab %>% data.frame() %>% dplyr::select(Df, `F`, P = Pr..F.)
  
  df %<>% rownames_to_column(var = 'Term')
  df$Term %<>% str_replace('Residuals', 'Residual')
  df$Term %<>% str_replace('Groups', test)
  df$`F` %<>% round(digits = 3)
  df$P %<>% round(digits = 3)
  
  total.df <- data.frame(Term = 'Total', Df = sum(df$Df))
  df %<>% bind_rows(total.df)
  
  df
}
permanova <- function(in.phy, method = 'bray', n.perm = 999, by = 'terms', form, block = 'none') {
  # This function accepts phyloseq objects as input, performs PERMANOVA,
  # and returns a data frame for visualization.
  # Complex formulae are accepted and relevant permutation options for vegan::adonis2 can be set.
  
  require(phyloseq)
  require(vegan)
  
  form <- as.formula(form)
  
  dist <- in.phy %>% phyloseq::distance(method = method)
  otu.tab <- in.phy@otu_table %>% data.frame()
  meta.tab <- in.phy@sam_data %>% data.frame()
  age <- meta.tab$age %>% unique()
  tree <- meta.tab$tree %>% unique()
  if(length(age) > 1){
    age <- 'All'
  }
  if(length(tree) > 1){
    tree <- 'All'
  }
  
  set.seed(666)
  perm <- how(nperm = n.perm)
  if(block == 'tree'){
    setBlocks(perm) <- with(meta.tab, tree) 
  }
  
  else if(block == 'group'){
    setBlocks(perm) <- with(meta.tab, group) 
  }
  
  else if(block == 'age'){
    setBlocks(perm) <- with(meta.tab, age) 
  }
  
  results <- adonis2(formula = form, data = meta.tab, by = by, permutations = perm)
  
  results %<>% data.frame() %>% rownames_to_column(var = 'Term') %>%
    dplyr::select(Term, Df, R2, `F`, P = Pr..F.) %>%
    mutate(`F` = round(`F`, 3),
           R2 = round(R2, 3),
           P = round(P, 3),
           .keep = 'unused')

  results$Age <- age
  results$Tree <- tree
  results
}
time.mantel <- function(age.1, age.2, group, spatial = F, dist = 'bray', method = 'spearman', n.perm = 999){
  # This function accepts phyloseq objects (representing different age classes and containing the same number of samples)
  # as input, performs a Mantel test, and returns a dataframe for visualization.
  # An exposure group option must be chosen, but relevant options can be passed to phyloseq::distance and vegan::mantel.
  # The arranging steps are crucial for ensuring that the matrices line up properly and that strata is properly set!
  # "tree" is hardcoded as a permutation stratum, so change the function if this isn't desired.
  
  require(phyloseq)
  
  sam.data.1 <- age.1 %>% sample_data() %>% data.frame() %>% arrange(tree.ht)
  sam.data.2 <- age.2 %>% sample_data() %>% data.frame() %>% arrange(tree.ht)
  # xyz <- sam.data.1 %>% select(x, y, z)
  
  # sample_data(age.1) <- sam.data.1 %>% sample_data()
  # sample_data(age.2) <- sam.data.2 %>% sample_data()
  
  # otu_table(age.1) <- otu_table(age.1) %>% data.frame() %>% .[order(row.names(sam.data.1)), ] %>% otu_table(taxa_are_rows = F)
  # otu_table(age.2) <- otu_table(age.2) %>% data.frame() %>% .[order(row.names(sam.data.2)), ] %>% otu_table(taxa_are_rows = F)
  
  otu.tab.1 <- otu_table(age.1) %>% data.frame() %>% .[match(rownames(sam.data.1), rownames(.)), ]
  otu.tab.2 <- otu_table(age.2) %>% data.frame() %>% .[match(rownames(sam.data.2), rownames(.)), ]
  
  young <- sam.data.1$age %>% unique()
  old <- sam.data.2$age %>% unique()
  ages <- paste0(young, ' to ', old)
  trees <- sam.data.1[sam.data.1$group == group, ]$tree

  samps.1 <- sam.data.1[sam.data.1$group == group, ] %>% rownames()
  samps.2 <- sam.data.2[sam.data.2$group == group, ] %>% rownames()

  # d1 <- distance(age.1, method = dist) %>% as.matrix()
  # d2 <- distance(age.2, method = dist) %>% as.matrix()
  
  d1 <- vegan::vegdist(otu.tab.1, method = dist) %>% as.matrix()
  d2 <- vegan::vegdist(otu.tab.2, method = dist) %>% as.matrix()
  
  # d3 <- dist(xyz) %>% as.matrix()

  d1 <- d1[rownames(d1) %in% samps.1,
           colnames(d1) %in% samps.1]
  # d1 <- d1[rownames(d1) %in% samps.1, ]
  # d1 <- d1[, colnames(d1) %in% samps.1]
  d2 <- d2[rownames(d2) %in% samps.2,
           colnames(d2) %in% samps.2]
  # d2 <- d2[rownames(d2) %in% samps.2, ]
  # d2 <- d2[, colnames(d2) %in% samps.2]
  # d3 <- d3[rownames(d3) %in% samps.1, colnames(d3) %in% samps.1]

  set.seed(666)
  results <- mantel(d1, d2, method = method, permutations = n.perm, strata = trees)
  # results <- mantel.partial(d1, d2, d3, method = method, permutations = n.perm, strata = trees)

  data.frame(Group = group, Transition = ages, p = results$signif, r = results$statistic)
  
}
sample.tree.hts <- function(sam.data.1, sam.data.2, otu.tab.1, otu.tab.2, group, dist, method){
  # Iterated function called by "boot.it.mantel". Randomly samples tree heights
  # shared by the input age classes to create simulated OTU tables for each.
  # Distance matrices are calculated from these tables, and distances involving samples
  # associated with the supplied exposure group are subset into final matrices.
  # The correlation between these matrices (each associated with an age class) is
  # calculated according the method supplied by the user and output to "boot.it.mantel".
  
  resamp <- sample(sam.data.1$tree.ht, replace = T)
  # resamp <- sample(sam.data.1[sam.data.1$group == group, ]$tree.ht, replace = T)
  resamp.2 <- resamp.1 <- data.frame(tree.ht = resamp)
  
  samps.1 <- sam.data.1 %>% filter(tree.ht %in% unique(resamp)) %>% .[.$group == group, ] %>% .$sample %>% sample(replace = T)
  samps.1 %<>% unique()
  samps.1.tree.ht <- sam.data.1 %>% filter(sample %in% samps.1) %>% .$tree.ht %>% unique()
  
  samps.2 <- sam.data.2 %>% filter(tree.ht %in% samps.1.tree.ht) %>% .[.$group == group, ] %>% .$sample
  samps.2 %<>% unique()
  
  resamp.1 %<>% left_join(otu.tab.1, by = 'tree.ht')
  resamp.1$tree.ht <- NULL
  resamp.1$sample <- paste0(resamp.1$sample, '.', 1:length(resamp.1$sample))
  resamp.1 %<>% column_to_rownames(var = 'sample')
  
  resamp.2 %<>% left_join(otu.tab.2, by = 'tree.ht')
  resamp.2$tree.ht <- NULL
  resamp.2$sample <- paste0(resamp.2$sample, '.', 1:length(resamp.2$sample))
  resamp.2 %<>% column_to_rownames(var = 'sample')
  
  d1 <- vegan::vegdist(resamp.1, method = dist) %>% as.matrix()
  d2 <- vegan::vegdist(resamp.2, method = dist) %>% as.matrix()
  
  d1 <- d1[str_remove(rownames(d1), '.[[:digit:]]+$') %in% samps.1,
           str_remove(colnames(d1), '.[[:digit:]]+$') %in% samps.1] %>% as.dist()
  # d1 <- d1[str_remove(rownames(d1), '.[[:digit:]]+$') %in% samps.1, ] %>% as.dist()
  # d1 <- d1[, str_remove(colnames(d1), '.[[:digit:]]+$') %in% samps.1] %>% as.dist()
  
  
  d2 <- d2[str_remove(rownames(d2), '.[[:digit:]]+$') %in% samps.2,
           str_remove(colnames(d2), '.[[:digit:]]+$') %in% samps.2] %>% as.dist()
  # d2 <- d2[str_remove(rownames(d2), '.[[:digit:]]+$') %in% samps.2, ] %>% as.dist()
  # d2 <- d2[ ,str_remove(colnames(d2), '.[[:digit:]]+$') %in% samps.2] %>% as.dist()
  
  cor(x = d1, y = d2, method = method)
}
boot.it.mantel <- function(age.1, age.2, group, dist = 'bray', method = 'spearman', n = 999, cores = 6){
  # This function accepts two phyloseq objects as input, each representing two different age classes.
  # Each object contains the same number of samples, and each sample in one object shares
  # a unique "tree.ht" value with only one other sample in the other object.
  # Calling "sample.tree.hts" as many times as desired by the user,
  # "boot.it.mantel" calculates the bootstrapped mean Mantel correlation and its associated confidence interval
  # for the exposure group supplied supplied by the user. This function can be parallelized as desired.
  
  # The output is a row of a data frame that can be bound with other rows and easily visualized.
  # The generation of means and confidence intervals for Mantel correlations allows correlations to be compared
  # across different transitions and exposure groups in a nonparametric fashion.
  
  require(parallel)
  
  sam.data.1 <- age.1 %>% sample_data() %>% data.frame() %>% arrange(tree.ht)
  sam.data.2 <- age.2 %>% sample_data() %>% data.frame() %>% arrange(tree.ht)
  
  otu.tab.1 <- age.1 %>% otu_table() %>% data.frame() %>% .[match(rownames(sam.data.1), rownames(.)), ] %>% rownames_to_column(var = 'sample')
  otu.tab.1$tree.ht <- sam.data.1$tree.ht
  
  otu.tab.2 <- age.2 %>% otu_table() %>% data.frame() %>% .[match(rownames(sam.data.2), rownames(.)), ] %>% rownames_to_column(var = 'sample')
  otu.tab.2$tree.ht <- sam.data.2$tree.ht
  
  sam.data.1 %<>% rownames_to_column(var = 'sample')
  sam.data.2 %<>% rownames_to_column(var = 'sample')
  
  young <- sam.data.1$age %>% unique()
  old <- sam.data.2$age %>% unique()
  ages <- paste0(young, ' to ', old)

  clust <- makeCluster(cores, setup_strategy = 'sequential')
  clusterEvalQ(clust, library(MASS))
  clusterEvalQ(clust, library(tidyverse))
  clusterEvalQ(clust, library(magrittr))
  clusterEvalQ(clust, source('code/functions.R'))
  clusterExport(clust, c('sam.data.1', 'sam.data.2',
                         'otu.tab.1', 'otu.tab.2'), envir = environment())
  clusterSetRNGStream(clust)
  boots <- parSapply(cl = clust, X = 1:n,
                     FUN = function(i){sample.tree.hts(sam.data.1, sam.data.2,
                                                       otu.tab.1, otu.tab.2,
                                                       group, dist, method)})
  stopCluster(clust)

  data.frame(Transition = ages,
             Group = group,
             mean = mean(boots),
             lci = coxed::bca(boots)[1],
             uci = coxed::bca(boots)[2])
  
  # list(sam.data.1 = sam.data.1, otu.tab.1 = otu.tab.1, sam.data.2 = sam.data.2, otu.tab.2 = otu.tab.2)
}

# Taxon responses
indicate <- function(phy.in, test = 'group', prop = 0.1, iterate = 999, max.p = 0.05){
  # This function accepts a phyloseq object, two different categorical variables (test = "group" or "age"),
  # a proportion-based filtering threshold (prop), an iteration requirement, and a minimum accepted p-value.
  # Indicator species analysis is performed
  
  require(labdsv)
  set.seed(666)
  
  phy.in %<>% filt.n.prune(prop = prop, depth = 0) %>%
    transform_sample_counts(fun = rel)
  
  otu.tab <- phy.in@otu_table %>% data.frame()
  meta.tab <- phy.in@sam_data %>% data.frame()
  tax.tab <- phy.in@tax_table %>% data.frame() %>% dplyr::select(Taxon = Genus_species) %>%
    rownames_to_column('OTU')
  tax.tab$Taxon %<>% str_replace('Rasutoria', 'Zasmidium')
  age <- meta.tab$age %>% unique()
  if(length(age) > 1){
    age <- 'All'
  }
  
  if(test == 'group'){
    ind <- indval(otu.tab, meta.tab$group, numitr = iterate)
    ind.val <- data.frame(Age = age,
                          Group = ind$maxcls,
                          IV = ind$indcls,
                          P = ind$pval %>% p.adjust('fdr'))
    ind.val$Group %<>% str_replace('1', 'Closed') %>% str_replace('2', 'Open')
    ind.val$Group %<>% as.factor()
  }
  
  else if(test == 'age'){
    ind <- indval(otu.tab, meta.tab$age, numitr = iterate)
    ind.val <- data.frame(Age = ind$maxcls,
                          IV = ind$indcls,
                          P = ind$pval %>% p.adjust('fdr'))
    ind.val$Age %<>% str_replace('1', 'A1') %>%
      str_replace('2', 'A2') %>%
      str_replace('3', 'A3') %>%
      str_replace('4', 'A4')
  }
  
  ind.val$Age %<>% as.factor()
  ind.val %<>% filter(P < max.p)
  ind.val %<>% rownames_to_column(var = 'OTU') %<>% left_join(tax.tab, by = 'OTU')
  ind.val$OTU %<>% str_replace('otu', 'OTU')
  ind.val %>% dplyr::select(OTU, Taxon, everything())
}
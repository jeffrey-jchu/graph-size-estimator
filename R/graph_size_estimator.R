#' Graph and Network Sampler
#'
#' A short function to perform simple random walk sampling on graph or network data.
#' Requires the igraph package.
#' @param nodes Graph or network data, in the form of a table of two columns representing edges between pairs of nodes: from (column 1) and to (column 2).
#' @param steps Number of steps to take before the next node is sampled in the random walk.
#' @param subset The number of nodes to be sampled in the random walk.
#' @details The output of the random walk can be accessed via two commands:
#' \itemize{
#' \item \code{graph.Random.Walk} returns a vector of nodes sampled;
#' \item \code{rw.Node.Degrees} returns a vector of the corresponding node degrees.
#' }
#' @keywords random walk
#' @keywords graph
#' @keywords network
#' @keywords sampling
#' @export
#' @examples
#' # Performs a random walk on the graph or network data,
#' # sampling each consectuive node until 1000 nodes have been sampled:
#' walk(nodes, 1, 1000)
#'
#' # Performs a random walk on the graph or network data,
#' # sampling nodes 3 steps apart until 5000 nodes have been sampled:
#' walk(nodes, 3, 5000)

rwalk <- function(nodes, steps, subset){
  # Read list of connections (edges) from 'nodes'
  # Generate a list of all unique nodes in the graph
  nodes.Unique <- unique(nodes$V1)

  # Use the igraph package to generate a random walk sequence on the graph, and node and degree lists
  # Generate an adjacency matrix representation of the graph
  node.Graph <- graph_from_data_frame(nodes, directed = FALSE, vertices = NULL)

  # Perform a random walk on the graph
  # Set the start node as a random node from the full graph
  node.Start <- sample((1:length(nodes.Unique)), 1)

  # Set the length of the random walk
  # (random walk length can be greater than the number of nodes sampled)
  # Use 'steps' and 'subset' to generate random walk length
  rw.Length = subset*steps
  # Generate the random walk, which outputs a vector of all the visited nodes
  graph.Random.Walk <- random_walk(node.Graph, node.Start, rw.Length, mode = 'all', stuck = 'return')
  # Obtain the degrees of each sampled node visited
  rw.Node.Degrees <- degree(node.Graph, v = graph.Random.Walk)

  # Convert the data frames of nodes and degrees into vectors
  graph.Random.Walk <- as.vector(graph.Random.Walk)
  rw.Node.Degrees <- as.vector(rw.Node.Degrees)

  # Generate a step sequence to obtain final random walk with steps
  # and gives vectors of only the nodes sampled and their corresponding degree values
  step.Sequence <- seq(1, rw.Length, steps)
  graph.Random.Walk <<- graph.Random.Walk[step.Sequence]
  rw.Node.Degrees <<- rw.Node.Degrees[step.Sequence]
}

#' Graph and Network Size Estimator
#'
#' A short function to compute estimates for the size of a graph or network, using small samples of data from the graph or network.
#' Requires the ggplot2 package.
#' @param node.File Vector of nodes sampled from the full graph or network.
#' @param degree.File Vector of the corresponding degrees of the nodes in \code{node.File}.
#' @param rw.Step Number of steps to take before the next node is sampled in the simple random walk.
#' @param size Sample size increment.
#' @param N True population value.
#' @param max.Rs.Division Maximum random sample divisor.
#' @param g Srivastava and Bhatnagar's constant, \code{g = 1} gives \code{N_B = N_A}.
#' @param wn.Constant Wither's and Nadarajah's constant (A).
#' @details Can be used in conjunction with \code{walk}, where \code{node.File = graph.Random.Walk},
#' and \code{degree.File = rw.Node.Degrees}.
#'
#' \code{size} can be computed by \deqn{size = (length(node.File)/rw.Step)/(max.Rs.Division + 1)}
#' @keywords graph
#' @keywords network
#' @keywords size
#' @keywords estimator
#' @export
#' @examples
#' # e.g. node.File and degree.File contain
#' # 2,000,000 nodes and degree values, respectively.
#' # The following computes estimates of
#' # the population size, and bias of the estimators, at sample sizes
#' # from n = 400 to 3600, in stpes of 400, with g = 10, and A = 0.5.
#' size.Est(node.File, degree.File, 500, 400, 42000000, 9, 10, 0.5)

size.Est <- function(node.File, degree.File, rw.Step, size, N, max.Rs.Division, g, wn.Constant){
  # Read nodes from the node data file
  node.List <- node.File
  # Read degrees from the degree data file
  degree.List <- degree.File
  # Random walk step length = rw.step
  # Sample size increments = size
  # True population size = N
  # Maximum random sample division limit = max.Rs.Division
  # Compute number of nodes in the data
  node.Length <- length(node.List)
  # Generate a sequence from 1 to node.Length
  node.Length.Sequence <- 1:(node.Length)

  # Initialise input parameters
  # Srivastava scalar = g
  # Withers-Nadarajah constant = wn.Constant
  # Reset Withers-Nadarajah summation
  total <- 0

  # Initialise outputs
  # Generate a matrix with columns equal to the division limit
  # and rows equal to the random walk step length - 1
  box.Gamma <- matrix(0, ncol = max.Rs.Division, nrow = (rw.Step - 1))
  # Generate matrices to store computations for each sample and estimator
  box.C <- box.Gamma # Collisions
  box.n <- box.Gamma # Sample sizes
  box.U <- box.Gamma # Unique nodes
  binomial.Coefficient <- box.Gamma # Binomial coefficient
  box.Dup <- box.Gamma # Duplicates
  box.Dmean <- box.Gamma #
  box.N <- box.Gamma # est N
  box.Nl <- box.Gamma # est Nl (Lu & Li)
  box.Na <- box.Gamma # est Na (Zellner)
  box.Nb <- box.Gamma # est Nb (Srivastava)
  box.Nc <- box.Gamma # est Nc (Voinov)
  box.Nd <- box.Gamma # est Nd (Withers & Nadarajah)

  # Random walk and collision counting loop
  # Loop from 1 to random walk step length - 1
  for (j in 1:(rw.Step - 1)) {
    # Generate sequence from j to number of nodes in data, in steps of rw.Step
    sample.Rw <- seq(j, node.Length, rw.Step)
    # Obtain nodes/degrees from data at 'sample.Rw' points in the data as a random walk
    sample.Nodes <- node.List[sample.Rw]
    sample.Degrees <- degree.List[sample.Rw]
    # Generate sample size from given parameters
    sample.Size <- floor((length(sample.Nodes))/(max.Rs.Division + 1))
    # Loop from 1 to maximum division limit
    for (random.Sample.Division in 1:(max.Rs.Division)) {
      # Reset random sample nodes and degrees
      random.Sample.Nodes <- c(0)
      random.Sample.Degrees <- c(0)
      # Assign random sample nodes and degrees
      random.Sample.Nodes <- sample.Nodes[(0*random.Sample.Division + 1):(sample.Size*(random.Sample.Division + 1))]
      random.Sample.Degrees <- sample.Degrees[(0*random.Sample.Division + 1):(sample.Size*(random.Sample.Division + 1))]
      n <- length(random.Sample.Nodes) # Obtain random sample size n
      tau <- sum(random.Sample.Degrees) # Compute sum of sample degrees
      # Compute sample mean of degrees of random walk sample
      # (estimates the asymptotic mean)
      d.Sample.Mean <- tau/n
      # Compute harmonic mean of sample degrees
      # (estimates the population mean of sample degrees)
      inverse.Sample.Degrees <- 1/(random.Sample.Degrees)
      d.Harmonic.Mean <- n/sum(inverse.Sample.Degrees)
      # Compute gamma^2 + 1 parameter from the sample
      gamma.Squared <- d.Sample.Mean/d.Harmonic.Mean
      # Assign random sample of nodes to a vector A
      A <- random.Sample.Nodes
      A <- as.matrix(A) # Convert vector of nodes to column matrix
      random.Nodes.Unique <- unique(A) # Find unique nodes in the sample
      unique.Length <- length(random.Nodes.Unique) # Obtain number of unique nodes
      a.Length <- length(A) # Obtain size of random sample of nodes
      duplicate.Length <- a.Length - unique.Length # Obtain number of duplicate nodes
      random.Nodes.Unique <- sort(random.Nodes.Unique) # Sort and order unique nodes
      # Collision and unique element counting
      # Vector of size equal to random sample size
      ff.Visits <- rep(0, a.Length)
      # Vector counting node 'i' visited 'x' number of times
      count.Visits <- match(A, random.Nodes.Unique) # Match uniques to random sample
      count.Visits <- sort(count.Visits) # Sort counts in ascending order
      # Generate count number from 1 to size of the random sample
      count.Number <- rep(0, a.Length)
      # Loop from 1 to size of the random sample
      for (i in 1:a.Length) {
        # Cycle through node 1, 2, ..., n
        count.Number[i] = count.Visits[i]
        # Count number of times node 'i' has been visited
        ff.Visits[count.Number[i]] = ff.Visits[count.Number[i]] + 1
      }
      # Vector counting number of nodes visited 'y' number of times, equal to number of unique nodes
      f <- rep(0, unique.Length)
      # Loop from 1 to size of random sample
      for (i in 1:a.Length) {
        # Check if node 'i' visted at least once
        if (ff.Visits[i] > 0)
          # Add to total of nodes visited 'y' number of times
          f[ff.Visits[i]] = f[ff.Visits[i]] + 1
      }
      # Remove entries where 0 nodes visited 'y' number of times
      # f = f[f != 0]
      # Null vector to reset and store collisions
      collision <- c(0)
      # Loop from 1 to length of max number of times a node was visited
      for (i in 1:length(f)) {
        # C = sum_{i}^{infty} \binom(n,2)*f
        collision = collision + (f[i]*(1/2)*i*(i - 1))
      }
      # Generate Katzir variation estimator Nk
      hat.Nk <- gamma.Squared*((n*(n - 1))/2)/(collision)
      # Assign (gamma^2 + 1); C; n; unique elements to matrices
      box.Gamma[j, random.Sample.Division] <- gamma.Squared
      box.C[j, random.Sample.Division] <- collision
      box.n[j, random.Sample.Division] <- n
      box.U[j, random.Sample.Division] <- unique.Length
      # Compute binomial coefficient using sample size n
      binomial.Coefficient[j, random.Sample.Division] <- n*(n - 1)/2

      # Assign number of duplicates; harmonic mean of degrees, and generate estimates
      # from the estimators (N; Nl) and assign to respective matrices
      box.Dup[j, random.Sample.Division] <- duplicate.Length
      box.Dmean[j, random.Sample.Division] <- d.Harmonic.Mean
      # box.Nk[j, random.Sample.Division] <- hat.Nk
      box.N[j, random.Sample.Division] <- binomial.Coefficient[j, random.Sample.Division]*gamma.Squared/(collision)
      box.Nl[j, random.Sample.Division] <- binomial.Coefficient[j, random.Sample.Division]*gamma.Squared/(collision + 1)
    }
  }

  # Compute additional parameters for Zellner estimator Na
  mean.Binomial.Coefficient <- apply(binomial.Coefficient, 2, mean)
  c.Bar <- apply(box.C, 2, mean)
  c.Bar2 <- c.Bar^2
  s.Na <- apply(box.C, 2, var)
  gamma.Squared.Na <- apply(box.Gamma, 2, mean)

  # Compute additional parameters for Srivastava estimator Nb
  # Set variance value the same as in Na estimator, as it is identical
  s.Nb <- s.Na
  # Gamma^2 is a universal value
  gamma.Squared.Nb <- gamma.Squared.Na
  # No additional parameters needed as denominator is c.Bar2 + g*(s.Nb/(rw.Step - 1))

  # Compute additional parameter for Voinov estimator Nc
  # Set variance value the same as in Na estimator, as it is identical
  s.Nc <- s.Na
  # Generate standard deviation from variance values
  root.S.Nc <- sqrt(s.Nc)
  # Gamma^2 is a universal value
  gamma.Squared.Nc <- gamma.Squared.Na
  # No additional parameters needed

  # Compute bias of estimators
  # Original estimator N
  mean.N <- apply(box.N, 2, mean)
  bias.N <- (mean.N - N)/N
  # Lu & Li estimator Nl
  mean.Nl <- apply(box.Nl, 2, mean)
  bias.Nl <- (mean.Nl - N)/N
  # Zellner estimator Na
  mean.Na <- gamma.Squared.Na*mean.Binomial.Coefficient*c.Bar/(c.Bar2 + (s.Na/(rw.Step - 1)))
  bias.Na <- (mean.Na - N)/N
  # Srivastava estimator Nb
  mean.Nb <- gamma.Squared.Nb*mean.Binomial.Coefficient*c.Bar/(c.Bar2 + g*(s.Nb/(rw.Step - 1)))
  bias.Nb <- (mean.Nb - N)/N
  # Voinov estimator Nc
  mean.Nc <- log(gamma.Squared.Nc)+
    log(mean.Binomial.Coefficient)+
    (log(sqrt(2*(rw.Step - 1)*pi)/root.S.Nc)+
       (((rw.Step - 1)*c.Bar2)/(2*s.Nc)))+
    pnorm(-(sqrt(rw.Step - 1)*c.Bar)/root.S.Nc, log.p = TRUE)
  mean.Nc <- exp(mean.Nc)
  bias.Nc <- (mean.Nc - N)/N

  # Compute additional parameters for Withers-Nadarajah estimator Nd
  s.Nd <- s.Na
  gamma.Squared.Nd <- gamma.Squared.Nc
  # Summation for Nd
  z <- floor((rw.Step^wn.Constant)) # Summation limit
  summation_total <- rep(0, max.Rs.Division)
  mean.Nd <- c(0)
  # For loop to compute summation part of Nd
  for (i in 1:max.Rs.Division) {
    summation <- c(0)
    for (k in 1:z) {
      d.Factorial <- choose(((rw.Step - 1)/2) + (k - 1) - 1, (k - 1))*(factorial(k - 1))
      summation[k] = (factorial(2*(k - 1)))*
        (1/(d.Factorial))*
        ((-(rw.Step - 1)*s.Nd[i]*(1/(4*rw.Step*c.Bar2[i])))^(k - 1))*
        ((1/factorial(k - 1)))
    }
    summation <- summation[is.finite(summation)]
    summation_total[i] <- sum(summation)
    mean.Nd[i] <- gamma.Squared.Nd[i]*mean.Binomial.Coefficient[i]*(1/c.Bar[i])*summation_total[i]
    bias.Nd <- (mean.Nd - N)/N
  }

  # Set plot parameters
  xx <- seq(size,(size*max.Rs.Division),size)
  yrange <- c(-0.5, 0.5) # Y-axis range
  xrange <- range(xx) # X-axis range
  low.Lim <- min(bias.N[is.finite(bias.N)],bias.Na[is.finite(bias.Na)],
                 bias.Nb[is.finite(bias.Nb)],bias.Nc[is.finite(bias.Nc)],
                 bias.Nd[is.finite(bias.Nd)],bias.Nl[is.finite(bias.Nl)])
  high.Lim <- max(bias.N[is.finite(bias.N)],bias.Na[is.finite(bias.Na)],
                  bias.Nb[is.finite(bias.Nb)],bias.Nc[is.finite(bias.Nc)],
                  bias.Nd[is.finite(bias.Nd)],bias.Nl[is.finite(bias.Nl)])

  # New plots
  # Create new data frame with relative bias of all estimators
  df <- data.frame(sample.Size = xx, rb = bias.N, rb.L = bias.Nl,
                   rb.A = bias.Na, rb.B = bias.Nb, rb.C = bias.Nc, rb.D = bias.Nd)
  df[df == "Inf"] <- NA # Replace infinite values with NA
  head(df) # Check header of dataframe

  # ggplot2 for relative bias vs log sample size (n/10)
  bias.Plot <<- ggplot(df, aes(x = sample.Size)) +
    geom_line(aes(y = rb, colour = "N (Original)")) +
    geom_line(aes(y = rb.L, colour = "N_L (Lu & Li)")) +
    geom_line(aes(y = rb.A, colour = "N_A (Zellner)")) +
    geom_line(aes(y = rb.B, colour = "N_B (Srvisatava)")) +
    geom_line(aes(y = rb.C, colour = "N_C (Voinov)")) +
    geom_line(aes(y = rb.D, colour = "N_D (Withers & Nadarajah)")) +
    scale_x_continuous(breaks = c(xx)) +
    labs(colour = "Estimator") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab(label = "Relative Bias") +
    xlab(label = "Sample Size n") +
    coord_cartesian(ylim = c(low.Lim, high.Lim)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")

  # Generate table of bias for all estimators
  estimators <- c(bias.N, bias.Nl, bias.Na, bias.Nb, bias.Nc, bias.Nd)
  estimators.Length <- length(estimators)/max.Rs.Division
  bias.Table <<- matrix(estimators, nrow = max.Rs.Division, ncol = estimators.Length)
  dimnames(bias.Table) <<- list(c(xx),c("N", "N_L", "N_A", "N_B", "N_C", "N_D"))

  # Generate table of estimates for all estimators
  estimates <- c(mean.N, mean.Nl, mean.Na, mean.Nb, mean.Nc, mean.Nd)
  estimates.Length <- length(estimates)/max.Rs.Division
  est.Table <<- matrix(estimates, nrow = max.Rs.Division, ncol = estimates.Length)
  dimnames(est.Table) <<- list(c(xx),c("N", "N_L", "N_A", "N_B", "N_C", "N_D"))
}

#
#
#
#
# 
# -- Damiano Fantini --
#
#

# --------------------------------------------------
# Module #1 - Classes and required methods
# --------------------------------------------------

#' An S4 class to represent a set of cells whose movements were tracked over time
#'
#' @slot images is a list of imported images
#' @slot proc_images is a list of processed images
#' @slot ops is a list keeping track of the operations executed on the object
#' @slot optimized is a list including results of the params auto-optimization (optional)
#' @slot centroids is a list of detected centroids
#' @slot positions is a data.frame of cell positions across stacks
#' @slot tracks is a numeric matrix of cell tracks
#' @slot params is a list of parameters used for the analysis
#' @slot stats is a list of stats computed for the cell tracks
#' @slot metadata is a list including labels about the image, and the experiment
#' 
#' @exportClass 
trackedCells <- setClass(
  #Name of the class
  "trackedCells",
  
  #define slots
  slots = list(
    images="list",
    proc_images="list",
    ops = "list",
    optimized = "list",
    centroids = "list", 
    positions = "data.frame", 
    tracks = "matrix", 
    params = "list",
    stats = "list",
    metadata = "list"
  ),
  
  # Make a runction to check compatibility
  validity = function(object)
  {
    if(length(object@images) > 1){
      return(TRUE)
    } else {
      return("Malformed Data")
    }
  }
)


setMethod("initialize", "trackedCells",
          function(.Object, x) {
            .Object <- callNextMethod(.Object)
            
            # check data and args
            chk <- (sum(names(x) %in% c("images", "dim", "attributes")) == 3)
            
            if(!chk)
              stop("Malformed Data")
            
            # Initialize
            P0 <- data.frame(row = 1, col = 1, tau = 1)
            P0 <- P0[-1, ]
            
            T0 <- matrix(NA, ncol = 4, nrow = 0)
            
            O0 <- list(images = 1, 
                       optimized_params = 0, 
                       custom_params = 0,
                       track = 0, 
                       stats = 0
            )
            
            #Assign
            .Object@images <- x
            .Object@proc_images <- list()
            .Object@ops <- O0
            .Object@optimized <- list()
            .Object@centroids <- list() 
            .Object@positions <- P0 
            .Object@tracks <- T0
            .Object@params <- list()
            .Object@stats <- list()
            .Object@metadata <- list(tiff_file = NA, 
                                     experiment = NA, 
                                     condition = NA, 
                                     replicate = NA)
            
            # return 
            .Object
            
          })

setMethod("show", signature(object = "trackedCells"), 
          function(object) {
            LNS <- list(
              " ~~~ An S4 trackedCells object ~~~",
              "",
              paste0("      + Num of images: ", object@images$dim$NumberImages),
              paste0("      + Optimized Params: ", ifelse(object@ops$optimized_params == 1, "Yes", "No")),
              paste0("      + Run w/ Custom Params: ", ifelse(object@ops$custom_params == 1, "Yes", "No")),
              paste0("      + Cells Tracked: ", ifelse(object@ops$track == 1, "Yes", "No")),
              paste0("      + Stats Computed: ", ifelse(object@ops$stats == 1, "Yes", "No")),
              "")
            
            for(lni in LNS){
              cat(lni, sep = "\n")  
            }
          })  


# --------------------------------------------------
# Module #2 - Custom and Convenience f(x)
# --------------------------------------------------

#' Return the Next Odd Integer
#'
#' Returns the smallest odd number bigger than the number(s) provided as the argument
#'
#' @param x a vector of class numeric
#'
#' @return a vector of class integer
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' 
#' @examples
#' cellTracker:::next_odd(2:5)
#'
#' @keywords internal
next_odd <- function(x) {
  y <- base::floor(x) + 1
  y <- base::ifelse(y %% 2 == 0, y + 1, y)
  return(y)
}


#' Shift Array Circularly
#'
#' Circularly shift the elements in an array by a user-defined number of positions. 
#' This emulates the behavior of the corresponding Matlab Circhsift function.
#' 
#' @param x a character, numeric, or logical vector with at least n + 1 elements
#' @param n an integer corresponding to the number of positions for the shift
#'
#' @return a vector corresponding to x (same size, same class), whose elements have been shifted 
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @examples
#' cellTracker:::circshift(1:10, -2)
#'
#' @keywords internal
circshift <- function(x, n = 1) {
  
  n <- as.integer(n[1])
  nn <- abs(n)
  
  if(is.vector(x) || is.list(x)) {
    len <- length(x)
  } else if (is.data.frame(x) || is.matrix(x)) {
    len <- nrow(x)
  }
  
  if (len == 1)
    return(x)
  
  if (nn >= len)
    stop("Bad n!")
  
  nu_idx <- 1:len
  if(n > 0) {
    nu_idx <- c((length(nu_idx) - nn + 1):length(nu_idx), 1:(length(nu_idx) - nn)) 
  } else if (n < 0) {
    nu_idx <- c((nn + 1):length(nu_idx), 1:nn)
  }
  
  if(is.vector(x) || is.list(x)) {
    y <- x[nu_idx]
  } else if (is.data.frame(x) || is.matrix(x)) {
    y <- x[nu_idx,]
  }
  return(y)
}


#' Add Dimension to a Molten Data Frame
#'
#' Creates a new (molten) data matrix where all elements of y ar added to each row of x. 
#' Each row in x is recycled for each element in y. Elements in y are added as the first column in the returned matrix.
#' 
#' @param x a matrix or data.frame with at least 1 row and 1 column. 
#' @param y a vector with elements that will be added to x
#'
#' @return a matrix with an extra column as compared to x
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @examples
#' cellTracker:::add_dimension(x = cbind(1:4, 4:1), y = c(9, 7))
#'
#' @keywords internal
add_dimension <- function(x, y) {
  w0 <- lapply(1:nrow(x), function(j1) {
    tmp <- x[j1, ]
    w1 <- lapply(1:length(y), function(j2) {
      c(y[j2], tmp)  
    })
    do.call(rbind, w1)
  })
  out <- do.call(rbind, w0)    
  out <- as.matrix(out)
  rownames(out) <- NULL
  colnames(out) <- NULL
  return(out)
}


#' Make Hypercube
#'
#' Creates a Molten Hypercube with a user-defined number of dimensions. The values supplied by the user
#' are used to fill each dimension. All possible combination of values are included in the resulting hyper cube. 
#'
#' @param vals vector of values used to fill the hyper cube 
#' @param dims integer indicating the number of dimensions. The resulting molden data frame will have a number of 
#' columns equal to dims
#'
#' @return Matrix corresponding to a molten hyper cube. The number of columns is equal to dims; 
#' the number of rows is equal to length(vals) ^ dims
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @examples
#' cellTracker:::make_hypercube(1:3, 3)
#'
#' @keywords internal
make_hypercube <- function(vals, dims) {
  xi <- as.matrix(cbind(vals))
  yi <- vals
  if (dims > 1){
    for(i in 1:(dims-1)) {
      xi <- add_dimension(x = xi, y = yi)
    }
  } else {
    return(NULL)
  }
  return(xi)
}


#' Clean And Reformat a Numeric Matrix
#'
#' Convert any matrix-lie object to a numeric Matrix, and coerces all the elements to integer. 
#' Row names and column names are removed
#'
#' @param x matrix or data.frame including numeric data (or data that can be coerced to integer)
#' 
#' @return numeric matrix with all its elements coerced to integer 
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @examples
#' tmp <- data.frame(A = 1:4, B=c(3.1, 2.8, 3.3, 9.1), C = FALSE)
#' cellTracker:::matfix(tmp)
#'
#' @keywords internal
matfix <- function(x) {
  xx <- as.data.frame(x)
  for ( j in 1:ncol(xx)) {
    xx[, j] <- as.integer(xx[, j])
  }
  colnames(xx) <- NULL
  rownames(xx) <- NULL
  return(as.matrix(xx))
}


#' Linear Convolution of a Numeric Matrix
#'
#' Performs a linear convoltion of a Numeric Matrix, using a user-suplied linear kernel.
#' The convlution can be executed in a column-wise fashion by setting the col.wise argument to TRUE. 
#' Alternatively, the convolution is performed in a row-wise fashion
#'
#' @param x numeric matrix that will be used as input for the convoltion; 
#' this matrix typically corresponds to an image where signal (high values) indicates the
#' presence of a cell or a cell-like particle
#' @param krnl numeric vector corresponding to th kernel that will be used for the convolution. 
#' Briefly, the kernel includes the weights that will be used to compute a weighted sum at each 
#' position of the input numeric matrix
#' @param col.wise logical; shall the linear convolution be performed in a column-wise or 
#' row-wise fashion
#'
#' @return Linearly convoluted numeric matrix. The resulting matrix has the same 
#' dimensions of the inut matrix 
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @examples
#' graphics::par(mfrow = c(1, 2))
#' tmp <- sapply(1:12, function(i) { (6 + abs(i - 6)) * c(1:10, 10:1) })
#' cnv.tmp <- cellTracker:::linear_conv2(tmp, c(-3, 0, 3))
#' graphics::image(tmp); graphics::image(cnv.tmp)
#' @importFrom graphics par image
#'
#' @keywords internal
linear_conv2 <- function(x, krnl, col.wise = TRUE) 
{
  
  # Adjust based on col.wise
  if (col.wise) {
    xx <- t(x)
  } else {
    xx <- x
  }
  
  # Enlarge x based on kernel size
  ncl <- ncol(xx)
  tmp.i <- sapply(1:floor(length(krnl)/2), function(w) {xx[,1]})
  tmp.f <- sapply(1:floor(length(krnl)/2), function(w) {xx[,ncl]})
  X <- cbind(tmp.i, xx, tmp.f)
  
  # Proceed with convolution
  Y <- do.call(rbind, lapply(1:nrow(X), function(ri) {
    sapply(1:(ncol(X) - length(krnl) + 1), function(ci) {
      xcoord <- ci:(ci+length(krnl)-1)
      tmp <- X[ri, xcoord]
      as.numeric(rbind(tmp) %*% cbind(krnl))
    })  
  })) 
  
  if (col.wise) 
    Y <- t(Y)
  
  return(Y)
}






# --------------------------------------------------
# Module #3 - Visualization and plotting
# --------------------------------------------------

#' Visualize a matrix image
#'
#' Shows an image representation of a numeric matrix. Typically, this is a non-negative numeric matrix, 
#' where signal (high values) corresponds to the presence of cells, or cell-like particles
#'
#' @param img_mtx numeric matrix corresponding to a image
#' @param col character vector corresponding to a valid color palette 
#' @param ... additional arguments will be passed to graphics::image()
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#' 
#' @keywords cellTracker
#'
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics image
#'
#' @keywords internal
visualize_img <- function(img_mtx, col = NULL, ...) 
{
  
  if(is.list(img_mtx)) {
    img_mtx <- img_mtx[[1]]
  }
  
  if (is.null(col)) {
    col <- grDevices::colorRampPalette(c("white", "blue4"))(100)
  }
  
  if(!is.matrix(img_mtx))
    stop("The IMG is not a matrix!")
  
  m <- nrow(img_mtx)
  n <- ncol(img_mtx)
  xx <- t(img_mtx[m:1, ])
  graphics::image(xx, col = col, ...)
}




#' Visualize Cells in an Image Stack
#'
#' Visualize objects that were identified as cells in a given image stack
#'
#' @param tc_obj a trackedCells object
#' @param stack index of the image stack to use
#' @param pnt.cex cex of the points drawn around cells
#' @param txt.cex cex of the text used to annotate cells
#' @param offset offset value for the annotation
#' @param main string used for the plot title, can be NULL= NULL
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @export
visualize_stack_centroids <- function(tc_obj, stack = 1, 
                                      pnt.cex = 1.2, txt.cex = 0.9, 
                                      offset = 0.18, main = NULL) {
  
  
  b <- tc_obj@proc_images$images[[stack]]
  cnt <- tc_obj@centroids[[stack]]
  
  if(is.null(main)){
    main <- paste0("Stack num. ", stack)
  }
  
  visualize_img(img_mtx = b, las = 1, main = main)
  visualize_cntr(centroids = cnt, width_px = ncol(b), height_px = nrow(b), 
                 pnt.cex = pnt.cex, txt.cex = txt.cex, offset = offset) 
  
}

#' Visualize Centroids
#' 
#' Annotates centroids over an image
#' 
#' @param centroids centroid data.frame
#' @param width_px width of the image in pixels
#' @param height_px height of the image in pixels
#' @param pnt.cex cex of the point (circle) drawn around each cell
#' @param txt.cex cex of the text used to annotate the image
#' @param offset offset for the text annotations
#' @param col color of the points, e.g. "red2"
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#' 
#' @importFrom graphics points text
#'
#' @keywords internal
visualize_cntr <- function(centroids, width_px, height_px, pnt.cex = 1.2, 
                           txt.cex = 0.9, offset = 0.18, col = "red2")
{
  cnt <- centroids
  points(x = ((cnt$col - 1) / (width_px - 1)),
         y = 1-((cnt$row - 1) / (height_px - 1)), 
         cex = pnt.cex, col = col)
  
  points(x = ((cnt$col - 1) / (width_px - 1)),
         y = 1-((cnt$row - 1) / (height_px - 1)), 
         cex = 1.2, col = col)
  
  text(x = ((cnt$col - 1) / (width_px - 1)),
       y = 1-((cnt$row - 1) / (height_px - 1)), 
       labels = 1:nrow(cnt), font = 4, 
       cex = txt.cex, col = col, 
       pos = 4, offset = offset)
  
  #return()
}






# --------------------------------------------------
# Module #4 - Core FastTracks f(x)
# --------------------------------------------------



#' Import Image from TIFF
#'
#' Import a .tif stack containing fluorescently labeled point particles to be tracked
#'
#' @param tiff_file path to a TIFF file to be read in
#' @param experiment string, a label to describe the experiment (optional). Can be NULL
#' @param condition string, a label to describe the experimental condition (optional). Can be NULL 
#' @param replicate string, a label to identify the replicate (optional). Can be NULL
#'  
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#' 
#' @importFrom tiff readTIFF
#'
#' @export
load_tif <- function(tiff_file, experiment = NULL, condition = NULL, replicate = NULL) 
{
  myIMG <- tiff::readTIFF(source = tiff_file, 
                          native = FALSE, 
                          all = TRUE, 
                          info = TRUE,  
                          as.is = TRUE)
  
  if (!is.list(myIMG))
    myIMG <- list(myIMG)
  
  if (is.null(experiment)) {
    experiment <- NA  
  } else {
    experiment <- tryCatch(as.character(experiment[1]), error = function(e) NA)
  }
  
  if (is.null(replicate)) {
    replicate <- NA  
  } else {
    replicate <- tryCatch(as.character(replicate[1]), error = function(e) NA)
  }
  
  if (is.null(condition)) {
    condition <- NA  
  } else {
    condition <- tryCatch(as.character(condition[1]), error = function(e) NA)
  }
  
  # num of images
  NumberImages <- length(myIMG)
  
  # m means width... should be n of cols
  mImage <- ncol(myIMG[[1]])
  
  # n means height... should be n of rows
  nImage <- nrow(myIMG[[1]])
  
  # Get image INFO
  InfoImage <- try({lapply(myIMG, attributes)}, silent = TRUE)
  
  # Get image matrices
  FinalImage <- try({lapply(myIMG, function(x) {
    sapply(1:ncol(x), function(ii) {as.numeric(x[,ii])})
  })}, silent = TRUE)
  
  #return(list(images = FinalImage, 
  #            dim = list(NumberImages  = NumberImages , width_m = mImage, height_n = nImage),
  #            attributes = InfoImage))
  img_list <-   list(images = FinalImage,
                     dim = list(NumberImages  = NumberImages , width_m = mImage, height_n = nImage),
                     attributes = InfoImage)
  
  Y <- new(Class = "trackedCells", img_list)
  
  # Attach labels
  Y@metadata <- list(tiff_file = sub("^.*[/]([^/]+$)", "\\1", tiff_file),
                     experiment = experiment, 
                     condition = condition, 
                     replicate = replicate)
    
  return(Y)
}



#' Validate Centroids
#'
#' Validate parameters used to identify cells in a image stack. A figure containing 
#' current image frame with identified particles labeled with circles and numberical tags is generated
#'
#' @param stack stack of images to be evaluated
#' @param slice index of the frame within the stack to be evaluated
#' @param lobject integer, length in pixels somewhat larger than a typical object (cell)
#' @param threshold the minimum brightness of a pixel that might be local maxima. NOTE: 
#' Make it big and the code runs faster but you might miss some particles.  
#' Make it small and you'll get everything and it'll be slow.
#' @param pnt.cex cex of the circle drawn around each cell
#' @param txt.cex cex of the text used for annotating cells
#' @param offset offset used for annotating cells  
#'
#' @return data.frame of centroid positions
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @importFrom graphics box axis points text
#'
#' @keywords internal
centroid_validation <- function(stack, slice, lobject, threshold, pnt.cex = 1.2, txt.cex = 0.85, offset = 0.18)
{
  a <- stack$images[[slice]]
  b <- bpass(image_array = a, lnoise = 1, lobject = lobject, threshold = threshold)
  pk = pkfnd(im = b, th = threshold, sz = lobject+1)
  cnt = cntrd(im = b, mx = pk, sz = lobject+1)
  
  visualize_img(b, axes = FALSE)
  graphics::box()
  my_xax <- ncol(b)
  my_yax <- nrow(b)
  
  my_xax <- unique(c(seq(1, my_xax, by = 100), my_xax))
  my_yax <- unique(c(seq(1, my_yax, by = 100), my_yax))
  
  axis(side = 1,at = ((my_xax - 1) / max(my_xax)), labels = my_xax)
  axis(side = 2,at = 1-((my_yax - 1) / max(my_yax)), labels = my_yax, las = 1)
  
  
  points(x = ((cnt$col - 1) / (ncol(b) - 1)),
         y = 1-((cnt$row - 1) / (nrow(b) - 1)), 
         cex = pnt.cex, col = "red2")
  
  points(x = ((cnt$col - 1) / (ncol(b) - 1)),
         y = 1-((cnt$row - 1) / (nrow(b) - 1)), 
         cex = 1.2, col = "red2")
  
  text(x = ((cnt$col - 1) / (ncol(b) - 1)),
       y = 1-((cnt$row - 1) / (nrow(b) - 1)), 
       labels = 1:nrow(cnt), font = 4, 
       cex = txt.cex, col = "red2", 
       pos = 4, offset = offset)
  
  return(cnt)
}



#' Perform a bandpass by convolving with an appropriate kernel
#'
#' Implements a real-space bandpass filter that suppresses pixel noise and long-wavelength 
#' image variations while retaining information of a characteristic size.
#' First, a lowpassed image is produced by convolving the original with a gaussian.  
#' Next, a second lowpassed image is produced by convolving the original with a 
#' boxcar function. By subtracting the boxcar version from the gaussian version, 
#' we are using the boxcar version to perform a highpass.
#' This code 'bpass.pro' is copyright 1997, John C. Crocker and 
#' David G. Grier.  It should be considered 'freeware'- and may be
#' distributed freely in its original form when properly attributed.  
#' 
#' @param image_array Numeric matrix corresponding to the image to be filtered
#' @param lnoise Characteristic lengthscale of noise in pixels.
#' @param lobject Integer length in pixels somewhat larger than a typical object
#' @param threshold By default, after the convolution, any negative pixels are reset to 0.  
#' Threshold changes the threshhold for setting pixels to 0. Positive values may be useful for removing
#' stray noise or small particles.  
#' 
#' @return Numeric matrix corresponding to the filtered image
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @keywords internal
bpass <- function(image_array, lnoise, lobject = NULL, threshold)
{
  
  cstm_normalize <- function(x) { x/sum(x) }
  
  # Make kernel (linear)
  gaussian_kernel <- cstm_normalize(exp(-(seq(-2.5, 2.5, length.out = ((10 * lnoise) + 1))^2)))
  
  if (!is.null(lobject))  
    boxcar_kernel <- cstm_normalize(rep(1, times = (2 * lobject) + 1))
  
  
  gconv <- linear_conv2(t(image_array), gaussian_kernel)
  gconv <- linear_conv2(t(gconv), gaussian_kernel)
  
  #visualize_img(image_array, col = colorRampPalette(c("white", "red2"))(100))
  #visualize_img(gconv)
  
  if (!is.null(lobject)) {
    bconv <- linear_conv2(t(image_array), boxcar_kernel)
    bconv <- linear_conv2(t(bconv), boxcar_kernel)
    filtered <- gconv - bconv
  } else {
    filtered <- gconv
  }
  
  # Zero out the values on the edges to signal that they're not useful.     
  lzero <- max(lobject, ceiling(5*lnoise))
  
  filtered[1:(round(lzero)),] <- 0
  filtered[(nrow(filtered) - round(lzero) + 1):nrow(filtered),] <- 0
  
  filtered[, 1:(round(lzero))] <- 0
  filtered[, (ncol(filtered) - round(lzero) + 1):ncol(filtered)] <- 0
  
  # Zero all values below threshold
  filtered[filtered < threshold] <- 0
  
  return(filtered)
}



#' Calculates Centroids
#' 
#' Calculates the centroid of bright spots to sub-pixel accuracy. Inspired by Grier & 
#' Crocker's feature for IDL, but greatly simplified and optimized for matlab, and then 
#' further ported to R. CREATED: Eric R. Dufresne, Yale University, Feb 4 2005.
#' 
#' @param im numeric matrix corresponding to the image to process
#' @param mx location of local maxima to pixel-levels accuracy
#' @param sz diamter of the window over which to average to calculate the centroid. should be big enough.
#' @param interactive numeric; if set to 1 (or any positive number), an image showing the 
#' computed centroids will be visualized
#'
#' @return a data.frame with 4 columns, containing, x, y, brightness, and 
#' the square of the radius of gyration for each cell.
#' 
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @importFrom graphics box axis title 
#'
#' @keywords internal
cntrd <- function(im, mx, sz, interactive = NULL) 
{
  # check interactive
  if(is.null(interactive))
    interactive <- 0
  
  # check sz
  if ((sz/2) == (floor(sz/2))) {
    sz <- sz + 1
    message("sz must be odd, like bpass")
    message(paste0("sz set to ", sz))
  }
  
  # check mx
  if (is.null(mx) || (!is.data.frame(mx)) || nrow(mx) == 0) {
    message('there were no positions inputted into cntrd. check your pkfnd theshold')
    return(NULL)
  }
  
  # Compute
  r <- (sz+1)/2
  
  # Create mask - window around trial location over which to calculate the centroid
  m <- 2*r
  x <- 0:(m-1) 
  cent <- (m-1)/2
  x2 <- (x-cent) ^ 2
  dst <- do.call(rbind, lapply(1:m, function(i){
    sqrt((i-1-cent)^2+x2)  
  }))
  
  ind <- dst < r
  
  msk <- sapply(1:ncol(ind), function(j) {as.numeric(ind[,j])})
  dst2 <- msk * (dst^2)
  ndst2 <- sum(dst2, na.rm = TRUE)
  
  nr <- nrow(im)
  nc <- ncol(im)
  
  # remove all potential locations within distance sz from edges of image
  ind <- mx$col > 1.5 * sz & mx$col < nc - 1.5*sz
  mx <- mx[ind, ]
  ind <- mx$row > (1.5*sz) & mx$row < nr - 1.5*sz
  mx <- mx[ind, ]
  
  nmx <- nrow(mx)
  
  # inside of the window, assign an x and y coordinate for each pixel
  xl <- do.call(rbind, lapply(1:(2*r), function(j) {
    (1:(2*r))  
  }))
  yl <- t(xl)
  
  #loop through all of the candidate positions
  pts <- list()
  for (i in 1:nmx) {
    #create a small working array around each candidate location, and apply the window function
    tmp <- msk * im[(mx$row[i] - floor(r) + 1):(mx$row[i] + floor(r)), 
                    (mx$col[i] - floor(r) + 1):(mx$col[i] + floor(r))]
    
    #calculate the total brightness
    norm <- sum(tmp, na.rm = TRUE)
    
    #calculate the weigthed average x location
    xavg <- sum(tmp * xl) / norm
    
    #calculate the weighted average y location
    yavg <- sum(tmp * yl) / norm
    
    #calculate the radius of gyration^2
    #rg=(sum(sum(tmp.*dst2))/ndst2);
    rg <- sum(tmp * dst2)/norm
    
    #concatenate it up
    pts[[length(pts) +1 ]] <- data.frame(row = mx$row[i]+yavg-r, 
                                         col = mx$col[i] + xavg - r, 
                                         norm = norm , 
                                         rg = rg)
    if (as.numeric(interactive) > 0) {
      visualize_img(img_mtx = tmp, axes = FALSE)
      box()
      axis(side = 1, at = seq(0, 1, length.out = ncol(tmp)), 
           labels = (mx$row[i] - floor(r) + 1):(mx$row[i] + floor(r)))
      axis(side = 2, at = seq(0, 1, length.out = nrow(tmp)), 
           labels = (mx$col[i] + floor(r)):(mx$col[i] - floor(r) + 1), las = 1)
      title(main = paste0("Cell number #", i), ylab = "y_pixel", 
            xlab = "x_pixel", font = 2, cex = 0.9)
      
      # Wait for user input from keyboard
      readline("Press Enter for Next Cell...")
    }
  }
  pts <- do.call(rbind, pts)
  rownames(pts) <- NULL
  return(pts)
}





#' Find Signal Peaks
#' 
#' Finds local maxima in an image to pixel level accuracy. This provides a rough guess 
#' of particle centers to be used by cntrd(). Inspired by the lmx subroutine of 
#' Grier and Crocker's. CREATED: Eric R. Dufresne, Yale University, Feb 4 2005.
#' 
#' @param im image to process, particle should be bright spots on dark background 
#' with little noise ofen an bandpass filtered brightfield image
#' @param th the minimum brightness of a pixel that might be local maxima. 
#' NOTE: Make it big and the code runs faster but you might miss some particles.
#' Make it small and you'll get everything and it'll be slow.
#' @param sz if your data is noisy, (e.g. a single particle has multiple local maxima), 
#' then set this optional keyword to a value slightly larger than the diameter of your blob. 
#' If multiple peaks are found withing a radius of sz/2 then the code will keep only the brightest.  
#' Also gets rid of all peaks within sz of boundary
#' 
#' @return a numeric data.frame with two columns, with the coordinates of local maxima 
#' 
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @importFrom reshape2 melt 
#' 
#' @keywords internal
pkfnd <- function(im, th, sz=NULL) 
{
  
  # # nested f(x)
  #my_melt <- function(data, varnames = NULL) {
  #  
  #  out.list <- list()
  #  for(ci in 1:ncol(data)) {
  #    for (ri in 1:nrow(data)) {
  #      out.list[[length(out.list) + 1]] <- 
  #        data.frame(Var1 = ri, Var2 = ci, value = data[ri, ci], 
  #                   stringsAsFactors = FALSE)
  #    }
  #  }
  #  OUT <- do.call(rbind, out.list)
  #  
  #  if (!is.null(varnames)) {
  #    colnames(OUT)[1:2] <- varnames
  #  }
  #  return(OUT)
  #}
  
  # find all the pixels above threshold
  ind <- im > th
  nr <- nrow(im) 
  nc <- ncol(im)
  
  # melt to have a list of points above threshold
  ind2 <- reshape2::melt(data = ind, varnames = c("row", "col"))
  ind2 <- ind2[ind2$value, ]
  
  # check each pixel above threshold to see if it's brighter than it's neighbors
  # THERE'S GOT TO BE A FASTER WAY OF DOING THIS.  I'M CHECKING SOME MULTIPLE TIMES,
  # BUT THIS DOESN'T SEEM THAT SLOW COMPARED TO THE OTHER ROUTINES, ANYWAY.
  keep <- list()
  for(i in 1:nrow(ind2)) {
    ri <- ind2$row[i]
    ci <- ind2$col[i]
    
    if (ri>1 & ri<nr & ci>1 & ci<nc) {
      z1 <- im[ri, ci]
      z2 <- as.numeric(im[(ri-1):(ri+1), (ci-1):(ci+1)])
      
      if (sum(z1 < z2, na.rm = TRUE) == 0) {
        keep[[length(keep) + 1]] <- i
      }
    }
  }
  
  # Next step
  npks <- length(keep)
  if(npks > 0) {
    keep <- do.call(c, keep)
    mx <- ind2[keep, ]
  } else {
    return(NULL)
  }
  
  # if size is specified, then get ride of pks within size of boundary (i.e., a margin from image edges)
  if (!is.null(sz) & npks>0) {
    # throw out all pks within sz of boundary;
    keep <- mx$row > sz & mx$row < (nr - sz + 1) & mx$col > sz & mx$col < (nc - sz + 1)
    mx<-mx[keep,]   
  }
  
  # prevent from finding peaks within size of each other
  npks <- nrow(mx)
  if (!is.null(sz) & npks > 1) {
    # CREATE AN IMAGE WITH ONLY PEAKS
    mask <- matrix(FALSE, nrow = nrow(im), ncol = ncol(im))
    for(i in 1:nrow(mx)) {
      mask[mx$row[i], mx$col[i]] <- TRUE
    }
    
    tmp <- matrix(0, nrow=nrow(im), ncol = ncol(im))
    tmp[mask] <- im[mask]
    
    # LOOK IN NEIGHBORHOOD AROUND EACH PEAK, PICK THE BRIGHTEST
    for (i in 1:nrow(mx)) {
      astep <- floor(sz/2)
      roi <- tmp[(mx$row[i] - astep):(mx$row[i] + astep), (mx$col[i] - astep):(mx$col[i] + astep)]
      imax <- which.max(roi)
      chkrow <- imax %% nrow(roi)
      myrow <- ifelse(chkrow == 0, nrow(roi), chkrow)
      mycol <- ifelse(chkrow == 0, floor(imax/nrow(roi)), floor(imax/nrow(roi)) + 1)
      mv <- roi[myrow, mycol]
      
      tmp[(mx$row[i] - astep):(mx$row[i] + astep), (mx$col[i] - astep):(mx$col[i] + astep)] <- 0
      tmp[(mx$row[i] - astep + myrow - 1), (mx$col[i] - astep + mycol - 1)] <- mv
    }
    
    ind <- tmp > th
    nr <- nrow(tmp) 
    nc <- ncol(tmp)
    
    # melt to have a list of points above threshold
    ind.f <- reshape2::melt(data = ind, varnames = c("row", "col"))
    ind.f <- ind.f[ind.f$value, 1:2]
    rownames(ind.f) <- NULL
    
    return(ind.f)
  } else {
    return(NULL)
  }
}


#' Build a Centroid Array
#' 
#' Create an array containing centroid data for particles identified in each frame
#' of the imported TIFF image stack
#' 
#' @param stack 3D matrix loaded to workspace from .tif stack
#' @param lobject Integer length in pixels somewhat larger than a typical object
#' @param threshold the minimum brightness of a pixel that might be local maxima
#'
#' @return data.frame of centroids, whith 4 columns corresponding to x-postion of centroid, 
#' y-postion of centroid, brightness, and square of the radius of gyration
#' 
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#' 
#'
#' @keywords internal
centroid_array <- function(stack, lobject, threshold)
{
  #determine the number of slices within the stack
  m <- stack$dim$width_m
  n <- stack$dim$height_n
  p <- stack$dim$NumberImages
  
  centroid <- list()
  for(i in 1:p) {
    a <-  stack$images[[i]]
    b <- bpass(image_array = a, 
               lnoise = 1, 
               lobject = lobject, 
               threshold = quantile(a, 0.25)) # maybe set to 0 or to threshold
    
    pk <- pkfnd(b, threshold, lobject+1)
    cnt <- cntrd(im = b, mx = pk, sz = lobject + 1)
    
    if(is.null(cnt) || nrow(cnt) < 1) {
      message(paste0('No centroids detectd in frame ', i, '...
                     \nCheck nuclei validation settings for this frame.'))  
    }
    centroid[[length(centroid) + 1]] <- cnt
    }
  return(centroid)
}


#' Detect Linear Paricle Diameters
#'
#' Estimates the diameters of particles in a numeric or logical vector
#'
#'
#' @param x numeric or logical vector
#'
#' @return data.frame including two columns: MPOS indicates the centroid position of a particle, 
#' and LEN indicates the diameter size
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @examples
#' cellTracker:::detect_radii(c(0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1))
#'
#' @keywords internal
detect_radii <- function(x) {
  
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  
  if(length(table(x)) > 2) {
    my.mean <- mean(x, na.rm = TRUE)
    x[x >= my.mean] <- 1
    x[x < my.mean] <- 0
  } 
  
  radii <- list()
  xx <- which(x == 1)
  
  if (length(xx) > 1) {
    LN <- 1
    p0 <- xx[1]
    p1 <- xx[1]
    
    for (j in 2:length(xx)) {
      if (xx[j] == (xx[(j-1)] + 1)) {
        LN <- LN + 1
        p1 <- xx[j]
        if (j == length(xx)) {
          yy <- data.frame(MPOS = mean(c(p0, p1)), LEN = LN)
          radii[[length(radii) + 1]]  <- yy
        }
      } else {
        yy <- data.frame(MPOS = mean(c(p0, p1)), LEN = LN)
        radii[[length(radii) + 1]]  <- yy
        p0 <- xx[j]
        p1 <- xx[j]
        LN <- 1
      }
    } 
    
  } else if (length(xx) == 1) {
    
    yy <- data.frame(MPOS = xx, LEN = 1)
    radii[[length(radii) + 1]]  <- yy
  } 
  
  
  if (length(radii) > 0 ) {
    radii <- do.call(rbind, radii)
  } else {
    radii <- NULL
  }
  
  return(radii)
}



#' Detect Paricle Diameters in a Numeric matrix
#'
#' Estimates the diameters of particles in a numeric matrix
#'
#'
#' @param x numeric matrix corresponding to a digital image
#' @param px.margin integer, number of pixels used as margin while searching/filtering for neighboring particles
#' @param quantile.val numeric, must be bigger than 0 and smaller than 1. 
#' Quantile for discriminating signal and background; only pixels with intensity higher than the corresponding 
#' quantile will count as signal while estimating particle diameters
#' @param plot logial, shall a histogram of the distribution of diameters be shown
#' 
#' @return list including summary stats and data about the particles found in the image
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @examples
#' a <- cbind(c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1), 
#'            c(1, 1, 0, 0, 0, 0, 0, 0, 1, 1), 
#'            c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#'            c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0), 
#'            c(0, 0, 0, 1, 1, 1, 0, 0, 0, 0))
#' graphics::image(a)
#' b <- estimate_diameter_range(a)
#' print(b$estim.cell.num)
#' print(b$raw)
#' 
#' @importFrom graphics image hist
#'
#' @export
estimate_diameter_range <- function(x, px.margin = 2, quantile.val = 0.99, plot = TRUE) {
  
  QNTS <- as.numeric(quantile(x, probs = quantile.val[1]))
  
  B <- x
  B[B < QNTS] <- 0
  B[B >= QNTS] <- 1
  
  rdds <- do.call(rbind, lapply(1:ncol(B), function(ii) {
    
    out <- detect_radii(B[,ii])
    if (!is.null(out)) {
      data.frame(RPOS = out$MPOS, CPOS = ii, LEN = out$LEN)
    }
  }))
  rdds$KEEP <- TRUE
  
  for (j in 1:nrow(rdds)) {
    if (rdds$KEEP[j]){
      
      tdm <- ( 2 * px.margin)  + rdds$LEN[j]
      ROWmin <- rdds$RPOS[j] - (0.5 * tdm)
      ROWmax <- rdds$RPOS[j] + (0.5 * tdm)
      COLmin <- rdds$CPOS[j] - (0.5 * tdm)
      COLmax <- rdds$CPOS[j] + (0.5 * tdm)
      
      keep <- rdds$RPOS >= ROWmin & rdds$RPOS <= ROWmax & 
        rdds$CPOS >= COLmin & rdds$CPOS <= COLmax & rdds$KEEP
      keep <- which(keep)
      keep <- keep[keep != j]
      if (length(keep) > 0) {
        curVal <- rdds$LEN[j]
        allValz <- rdds$LEN[keep]
        
        if (sum(curVal > allValz) == length(allValz)) {
          rdds$KEEP[keep] <- FALSE
        } else {
          rdds$KEEP[j] <- FALSE
        }
      }
    }
  }
  
  FINL <- rdds[rdds$KEEP,]
  
  yy <- list(estim.cell.num = sum(FINL$KEEP),
             q50.diam = median(FINL$LEN, na.rm = TRUE), 
             q75.diam = as.numeric(quantile(FINL$LEN, na.rm = TRUE, probs = 0.75)),
             q90.diam = as.numeric(quantile(FINL$LEN, na.rm = TRUE, probs = 0.90)),
             q95.diam = as.numeric(quantile(FINL$LEN, na.rm = TRUE, probs = 0.95)),
             raw = FINL)
  
  if (plot) {
    try(hist(FINL$LEN, breaks = seq(min(FINL$LEN, na.rm = TRUE), 
                                    max(FINL$LEN, na.rm = TRUE), length.out = 20),
             xlab = "Particle Diameter", las = 1, main = "Diam. Distribution", 
             col = "aquamarine3"), silent = TRUE); box()
  }
  
  return(yy)
}




#' Track cells
#'
#' Constructs n-dimensional trajectories from a scrambled list of particle 
#' coordinates determined at discrete times (e.g. in consecutive image frames)
#'
#'
#' @param xyzs an array listing the scrambled coordinates and data of the 
#' different particles at different times
#' @param maxdisp an estimate of the maximum distance that a particle would 
#' move in a single time interval
#' @param params a list containing a few tracking parameters that are needed 
#' for the analysis
#' 
#'
#' @return data.frame including cell tracks data
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @keywords internal
track <- function(xyzs, maxdisp, params) 
{
  #% ; maxdisp: an estimate of the maximum distance that a particle 
  #% ;     would move in a single time interval.(see Restrictions)
  #%  OPTIONAL INPUT:
  #%   param:  a structure containing a few tracking parameters that are
  #%       needed for many applications.  If param is not included in the
  #%       function call, then default values are used.  If you set one value
  #%       make sure you set them all:
  #% ;         param.mem: this is the number of time steps that a particle can be
  #% ;             'lost' and then recovered again.  If the particle reappears
  #% ;             after this number of frames has elapsed, it will be
  #% ;             tracked as a new particle. The default setting is zero.
  #% ;             this is useful if particles occasionally 'drop out' of
  #% ;             the data.
  #% ;         param.dim: if the user would like to unscramble non-coordinate data
  #% ;             for the particles (e.g. apparent radius of gyration for
  #% ;             the particle images), then positionlist should
  #% ;             contain the position data in positionlist(0:param.dim-1,*)
  #% ;             and the extra data in positionlist(param.dim:d-1,*). It is then
  #% ;             necessary to set dim equal to the dimensionality of the
  #% ;             coordinate data to so that the track knows to ignore the
  #% ;             non-coordinate data in the construction of the 
  #% ;             trajectories. The default value is two.
  #% ;         param.good: set this keyword to eliminate all trajectories with
  #% ;             fewer than param.good valid positions.  This is useful
  #% ;             for eliminating very short, mostly 'lost' trajectories
  #% ;             due to blinking 'noise' particles in the data stream.
  #%;          param.quiet: set this keyword to 1 if you don't want any text
  #% ; OUTPUTS:
  #  % ; result:  a list containing the original data rows sorted 
  #% ;     into a series of trajectories.  To the original input 
  #% ;     data structure there is appended an additional column 
  #% ;     containing a unique 'id number' for each identified 
  #% ;     particle trajectory.  The result array is sorted so 
  #% ;     rows with corresponding id numbers are in contiguous 
  #% ;     blocks, with the time variable a monotonically
  #% ;     increasing function inside each block.  For example:
  #  % ;     
  #% ;
  #% ;     NB: for t=1 in the example above, one particle temporarily
  #% ;     vanished.  As a result, the trajectory id=1 has one time
  #% ;     missing, i.e. particle loss can cause time gaps to occur 
  #% ;     in the corresponding trajectory list. In contrast:
  #  % ;
  #% ;     IDL> res = track(pos,5)
  #% ;
  #% ;     track will return the result 'res'
  #% ;         (x)      (y)      (t)          (id)
  #% ;     res = 15.1000      22.6000      0.00000      0.00000
  #% ;                   3.60000      5.00000      0.00000      1.00000
  #% ;               4.10000      5.50000      1.00000      1.00000
  #% ;               6.20000      4.30000      2.00000      1.00000
  #% ;               15.9000      20.7000      2.00000      2.00000
  #% ; 
  #% ;     where the reappeared 'particle' will be labelled as new
  #% ;     rather than as a continuation of an old particle since
  #% ;     mem=0.  It is up to the user to decide what setting of 
  #% ;     'mem' will yeild the highest fidelity .
  #% ; 
  #% ; SIDE EFFECTS:
  #  % ; Produces informational messages.  Can be memory intensive for
  #% ; extremely large data sets.
  #% ; RESTRICTIONS:
  #  % ; maxdisp should be set to a value somewhat less than the mean 
  #% ; spacing between the particles. As maxdisp approaches the mean
  #% ; spacing the runtime will increase significantly. The function 
  #% ; will produce an error message: "Excessive Combinatorics!" if
  #% ; the run time would be too long, and the user should respond 
  #% ; by re-executing the function with a smaller value of maxdisp.
  #% ; Obviously, if the particles being tracked are frequently moving
  #% ; as much as their mean separation in a single time step, this
  #% ; function will not return acceptable trajectories.
  #% ; PROCEDURE:
  #  % ; Given the positions for n particles at time t(i), and m possible
  #% ; new positions at time t(i+1), this function considers all possible 
  #% ; identifications of the n old positions with the m new positions,
  #% ; and chooses that identification which results in the minimal total
  #% ; squared displacement. Those identifications which don't associate
  #% ; a new position within maxdisp of an old position ( particle loss )
  #% ; penalize the total squared displacement by maxdisp^2. For non-
  #% ; interacting Brownian particles with the same diffusivity, this
  #% ; algorithm will produce the most probable set of identifications 
  #% ; ( provided maxdisp >> RMS displacement between frames ).
  #% ; In practice it works reasonably well for systems with oscillatory,
  #% ; ballistic, correlated and random hopping motion, so long as single 
  #% ; time step displacements are reasonably small.  NB: multidimensional
  #% ; functionality is intended to facilitate tracking when additional
  #% ; information regarding target identity is available (e.g. size or 
  #% ; color).  At present, this information should be rescaled by the
  #% ; user to have a comparable or smaller (measurement) variance than 
  #% ; the spatial displacements.
  #% ;
  
  # Initialize and stuff
  warn_log <- list()
  
  warn_message <- function(warn_log, quiet = FALSE) {
    
    warn_cycles <- NULL
    
    if (is.list(warn_log) && length(warn_log) > 0) {
      
      warn_cycles <- sort(unique(do.call(c, warn_log)))
      
      if (!quiet) {
        message(paste0("Difficult combinatorics encountered while processing slide(s): ", 
                       paste(warn_cycles, collapse = ", "), "."))
      }
    }
    return(warn_cycles)
  } 
  
  
  
  dd <- ncol(xyzs)
  
  # use default parameters if none given
  # if nargin==2
  # default values
  memory_b <- 0
  goodenough <- 0
  dim <- dd - 1
  quiet <- FALSE
  force_exec <- FALSE
  
  
  if(!is.null(params)) {
    if(is.list(params)) {
      
      if(!is.null(params$memory_b) && is.numeric(params$memory_b)) {
        memory_b <- params$memory_b
      }
      
      if(!is.null(params$goodenough) && is.numeric(params$goodenough)) {
        goodenough <- params$goodenough
      }
      
      if(!is.null(params$dim) && is.numeric(params$dim)) {
        dim <- params$dim
      }
      
      if(!is.null(params$quiet) && is.logical(params$quiet)) {
        quiet <- params$quiet
      }
      
      
      if(!is.null(params$force_exec) && is.logical(params$force_exec)) {
        force_exec <- params$force_exec
      }
      
    }
  }
  
  # % checking the input time vector
  # THis should be monotonically not-decreasing and not identical
  tau <- xyzs[, dd]
  st <- tau[2:length(tau)] - tau[1:(length(tau) - 1)]
  
  if (sum(st < 0) > 0) {
    message("", appendLF = TRUE)
    message("The time vector (tau) is not ordered")  
    return(NULL)
  }
  
  if (length(unique(tau)) == 1) {
    message("", appendLF = TRUE)
    message('All positions are at the same time... go back!')  
    return(NULL)
  }
  
  #--remove if useless
  info <- 1
  w <- which(st > 0)
  z <- length(w)
  z <- z + 1
  
  # % partitioning the data with unique times
  # the first two lines were skipped in the original file
  # they are included here for completeness
  #res = unq(t);
  # implanting unq directly
  
  indices <- which(tau - circshift(tau, -1) != 0)
  count <- length(indices)
  
  if (count > 0) {
    res <- indices  
  } else{
    res = length(tau)-1  
  }
  
  res <- c(1, res, length(tau))
  ngood <- res[2] - res[1] + 1
  eyes <- 1:ngood
  pos <- xyzs[eyes, 1:dim]
  istart <- 2
  n <- ngood;
  
  zspan <- 50;
  if (n > 200) {
    zspan <- 20  
  } 
  
  if (n > 500){
    zspan <- 10  
  } 
  
  # initialize a matrix with -1
  resx <- matrix((-1), nrow = zspan, ncol = n)
  
  # initialize a second matrix with -1
  bigresx <- matrix((-1), nrow = z, ncol = n)
  mem <- matrix(0, nrow = n, ncol = 1)
  #%  whos resx
  #%  whos bigresx
  uniqid <- 1:n;
  maxid <- n;
  
  # initialize olis
  #olist <- data.frame(x= 0, y = 0)
  #olist <- c(0,0)
  olist <- list()
  
  if (goodenough > 0) {
    dumphash <- matrix(0, nrow = n, ncol = 1)
    nvalid <- matrix(1, nrow = n, ncol = 1)
  } 
  
  #%  whos eyes;
  resx[1,] <- eyes
  
  #% setting up constants
  maxdisq <- maxdisp^2
  
  #% (Little) John calls this the setup for "fancy code" ???
  # Robin replies: Fancy? Where? You got to be kidding, man!!!
  
  notnsqrd <- (sqrt(n*ngood) > 200) && (dim < 7)
  
  if (notnsqrd) {
    #%;   construct the vertices of a 3x3x3... d-dimensional hypercube
    numbs <- 0:2
    cube <- make_hypercube(vals = numbs, dims = dim)
    
    #%   calculate a blocksize which may be greater than maxdisp, but which
    #%   keeps nblocks reasonably small.  
    volume <- 1
    for (d in 0:(dim-1)) {
      minn <- min(xyzs[w, (d+1)])
      maxx = max(xyzs[w, (d+1)])
      volume <- volume * (maxx-minn)
    }
    
    # volume;
    blocksize <- max( c(maxdisp,((volume)/(20*ngood))^(1.0/dim)) )
  }
  
  
  ### %   Start the main loop over the frames.
  for (i in istart:z){
    
    #message(paste0("i=", i))
    ispan <- ((i-1) %% zspan) + 1
    # %disp(ispan)
    # % get new particle positions
    m <- res[(i+1)] - res[(i)]
    # res[i]
    eyes <- 1:m
    eyes <- eyes + res[i]
    
    if (m > 0) {
      xyi <- xyzs[eyes, 1:dim]
      found <- matrix(0, nrow = m, ncol = 1)
      
      # % THE TRIVIAL BOND CODE BEGINS   
      if (notnsqrd) {
        
        # %Use the raster metric code to do trivial bonds
        
        #% construct "s", a one dimensional parameterization of the space 
        #% which consists of the d-dimensional raster scan of the volume.)
        
        abi <- matfix(xyi/blocksize)
        abpos <- matfix(pos/blocksize)
        si <- matrix(0, nrow = m, ncol = 1)
        spos <- matrix(0, nrow = n, ncol = 1)
        dimm <- matrix(0, nrow=dim, ncol=1)
        coff <- 1
        
        for (j in 1:dim){
          minn <- min(c(as.numeric(abi[,j]), 
                        as.numeric(abpos[, j])), na.rm = TRUE)
          maxx <- max(c(as.numeric(abi[,j]), 
                        as.numeric(abpos[,j])), na.rm = TRUE)
          abi[, j] <- abi[, j] - minn
          abpos[,j] <- abpos[,j] - minn
          dimm[j,1] <- maxx-minn + 1
          si <- si + abi[,j] * coff
          spos <- spos + abpos[,j]*coff
          coff <- dimm[j,1]*coff
        }
        nblocks <- coff
        #% trim down (intersect) the hypercube if its too big to fit in the
        #% particle volume. (i.e. if dimm(j) lt 3)
        
        cub <- cube
        deg <- which(dimm[,1] < 3)
        if (length(deg) > 0) {
          for (j in 0:(length(deg)-1)){
            cub <- cub[which(cub[, deg[j+1]] < dimm[deg[j+1],1]) ,]
          }
        } 
        
        # % calculate the "s" coordinates of hypercube (with a corner @ the origin)
        scube <- matrix(0, nrow = nrow(cub), ncol=1)
        coff <- 1
        for (j in 1:dim){
          scube <- scube + (cub[,j] * coff)
          coff <- coff*dimm[j, 1]      
        }
        
        # % shift the hypercube "s" coordinates to be centered around the origin
        coff <- 1
        for (j in 1:dim){
          if (dimm[j, 1] > 3) {
            scube <- scube - coff
          }
          coff <- dimm[j, 1] * coff
        }
        scube <- (scube + nblocks) %% nblocks
        
        # get the sorting for the particles by their "s" positions.
        ed <- sort(si)
        isort <- order(si)
        
        #% make a hash table which will allow us to know which new particles
        #% are at a given si.
        strt <- matrix((-1), nrow = nblocks, ncol = 1)
        fnsh <- matrix(0, nrow = nblocks, ncol = 1)
        h <- which(si == 0)
        lh <- length(h)
        if (lh > 0) {
          si[h] <- 1  
        }
        
        for (j in 1:m){
          if (strt[si[isort[j]], 1] == (-1)){
            strt[si[isort[j]],1] <- j
            fnsh[si[isort[j]], 1] <- j
          } else {
            fnsh[si[isort[j]], 1] <- j
          }
        }
        if (lh > 0) {
          si[h] <- 0   
        }
        
        
        coltot <- matrix(0, nrow = m, ncol = 1)
        rowtot <- matrix(0, nrow = n, ncol = 1)
        which1 <- matrix(0, nrow = n, ncol = 1)
        
        for (j in 1:n){
          
          map <- matfix(-1)
          
          scub_spos <- scube + spos[j];
          s <- scub_spos %% nblocks
          whzero <- which(s == 0 )
          if (length(whzero > 0)){
            nfk <- which(s !=0 )
            s <- s[nfk]
          }
          
          w <- which(strt[s, 1] != (-1))
          
          ngood <- length(w)
          ltmax <- 0
          if (ngood != 0){
            
            s <- s[w]
            for (k in 1:ngood){
              map = c(map, isort[strt[s[k]]:fnsh[s[k]]])
            }
            map <- map[2:length(map)]
            #%                     if length(map) == 2
            #%                         if (map(1) - map(2)) == 0
            #%                             map = unique(map);
            #%                          end
            #%                     end
            #%   map = map(umap);
            #%end
            #% find those trival bonds
            distq <- matrix(0, nrow=length(map), ncol=1)
            for (d in 1:dim){     
              distq <- distq + (xyi[map,d] - pos[j,d])^2
            }
            ltmax <- distq < maxdisq
            
            rowtot[j, 1] <- sum(ltmax)
            
            if (rowtot[j] >= 1){ 
              w <- which(ltmax == 1)
              coltot[map[w], 1] <- coltot[ map[w], 1] +1
              which1[j, 1] <- map[ w[1]]
            }
          }
        }
        
        
        ntrk <- matfix(n - sum(rowtot == 0))
        
        w <- which(rowtot == 1)
        ngood <- length(w)
        
        
        if (ngood != 0) { 
          ww <- which(coltot( which1[w] ) == 1);
          ngood <- length(ww)
          if (ngood != 0){ 
            # %disp(size(w(ww)))
            resx[ispan, w[ww]] <- eyes[which1[w[ww]]]
            found[which1[w[ww]]] <- 1
            rowtot[w[ww]] = 0;
            coltot[which1[w[ww]]] <- 0
          }
        }
        
        labely <- which(rowtot > 0)
        ngood <- length(labely)
        if (ngood != 0){ 
          labelx <- which(coltot > 0)
          
          nontrivial <- 1
        } else {
          nontrivial <- 0
        }
        
      } else { 
        
        # % THE TRIVIAL BOND {else} block CODE BEGINS   
        #%   or: Use simple N^2 time routine to calculate trivial bonds      
        
        #% let's try a nice, loopless way!
        #% don't bother tracking perm. lost guys.
        wh <- which(pos[,1] >= 0)
        ntrack <- length(wh)
        if (ntrack == 0){ 
          message('There are no valid particles to track!')
          break()
        }
        
        # yma initialization was added
        xmat <- matrix(0, nrow = ntrack, ncol = m)
        ymat <- matrix(0, nrow = ntrack, ncol = m)
        
        count <- 0
        for (kk in 1:ntrack) {
          for (ll in 1:m) {
            xmat[kk,ll] <- count
            count <- count+1
          }
        }
        count <- 0
        
        # if there are not enough cols or rows, add them and set to 0
        if (nrow(ymat) < m) {
          TMP <- matrix(0, nrow = (m - nrow(ymat)), ncol = ncol(ymat)) 
          ymat <- rbind(ymat, TMP)
        }
        if (ncol(ymat) < ntrack) {
          TMP <- matrix(0, nrow = nrow(ymat), ncol = (ntrack - ncol(ymat))) 
          ymat <- cbind(ymat, TMP)
        }
        
        for (kk in 1:m) {
          for (ll in 1:ntrack) {
            ymat[kk,ll] <- count
            count <- count+1
          }
        }
        
        xmat <- (xmat %% m) + 1
        ymat <- t((ymat %% ntrack) +1)
        lenxn <- nrow(xmat)
        lenxm <- ncol(xmat)
        #%            whos ymat
        #%            whos xmat
        #%            disp(m)
        
        for (d in 1:dim) {
          x <- xyi[,d]
          y <- pos[wh,d]
          
          xm <- sapply(1:ncol(xmat), function(jj) {
            tcljj <- xmat[, jj]
            x[tcljj]
          })
          
          #ym <- y[ymat[1:lenxn, 1:lenxm]]
          tmpymat <- ymat[1:lenxn, 1:lenxm]
          ym <- sapply(1:ncol(tmpymat), function(jj) {
            tcljj <- tmpymat[, jj]
            y[tcljj]
          })
          
          if (nrow(xm) != nrow(ym) || ncol(xm) != ncol(ym)) {
            xm <- t(xm)
          }
          
          if (d == 1) {
            dq <- (xm -ym)^2
            #%dq = (x(xmat)-y(ymat(1:lenxn,1:lenxm))).^2;
          } else {
            dq <- dq + (xm-ym)^2
            #%dq = dq + (x(xmat)-y(ymat(1:lenxn,1:lenxm)) ).^2;
          }
        }
        
        ltmax <- 1 * (dq < maxdisq)
        
        #% figure out which trivial bonds go with which
        
        rowtot <- matrix(0, nrow = n, ncol = 1)
        rowtot[wh, 1] <- apply(ltmax, 1, sum)
        
        if (ntrack > 1) { 
          coltot <- apply(ltmax, 2, sum, na.rm = TRUE)
        } else {
          coltot <- ltmax
        }
        which1 <- matrix(0, nrow = n, ncol = 1)
        for (j in 1:ntrack) { 
          mx  <- max(ltmax[j, ], na.rm = TRUE)
          w <- which.max(ltmax[j, ])
          which1[wh[j]] <- w
        }
        
        ntrk <- matfix( n - sum(rowtot == 0))
        w <- which( rowtot == 1)
        ngood <- length(w)
        if (ngood != 0) {
          ww <- which(coltot[which1[w]] == 1)
          ngood <- length(ww)
          if (ngood != 0) { 
            resx[ ispan, w[ww] ] <- eyes[ which1[w[ww]]]
            found[which1[ w[ww]]] <- 1
            rowtot[w[ww]] <- 0
            coltot[which1[w[ww]]] <- 0
          }
        }
        
        labely <- which(rowtot > 0)
        ngood <- length(labely)
        
        if (ngood != 0) {
          labelx <- which(coltot > 0)
          nontrivial <- 1
        } else {
          nontrivial <- 0
        }
      }
      
      # %THE TRIVIAL BOND CODE ENDS
      
      if (nontrivial == 1){
        
        xdim <- length(labelx)
        ydim <- length(labely)
        
        #%  make a list of the non-trivial bonds            
        
        bonds <- list()
        bondlen <- list()
        
        for (j in 1:ydim) {
          distq <- matrix(0, nrow = xdim, ncol = 1)
          
          for (d in 1:dim) {
            #%distq
            distq <- distq + cbind((xyi[labelx,d] - pos[labely[j],d])^2)
            #%distq    
          }
          
          w <- which(distq < maxdisq) - 1
          ngood <- length(w)
          newb <- rbind(w, rep(j, times = ngood))
          
          bonds[[(length(bonds) + 1)]] <- t(newb)
          bondlen[[(length(bondlen) + 1)]] = distq[w + 1]
        }
        
        bonds <- do.call(rbind, bonds)
        bondlen <- do.call(c, bondlen)
        
        numbonds <- length(bonds[,1])
        mbonds <- bonds;
        #max([xdim,ydim]);
        
        
        if (max(c(xdim,ydim)) < 4){
          nclust <- 1
          maxsz <- 0
          mxsz <- xdim
          mysz <- ydim
          bmap <- matrix((-1), nrow = length(bonds[,1]) + 1, 1)
          
        } else {
          
          #  %   THE SUBNETWORK CODE BEGINS            
          lista <- matrix(0, nrow = numbonds, ncol = 1)
          listb <- matrix(0, nrow = numbonds, ncol = 1)
          nclust <- 0
          maxsz <- 0
          thru <- xdim
          
          while (thru != 0) {
            #%  the following code extracts connected 
            #%   sub-networks of the non-trivial 
            #%   bonds.  NB: lista/b can have redundant entries due to 
            #%   multiple-connected subnetworks      
            
            w <- which(bonds[, 2] >= 0)
            
            lista[1] = bonds[w[1],2]
            listb[1] = bonds[w[1],1]
            bonds[w[1], ] <- (-1) * (nclust+1)
            # bonds;
            adda <- 1
            addb <- 1
            donea <- 0
            doneb <- 0
            if ((donea != adda) || (doneb != addb)){
              true <- FALSE
            } else {
              true <- TRUE   
            }
            
            while (!true){
              
              if (donea != adda) {
                w <- which(bonds[,2] == lista[donea+1])
                ngood <- length(w)
                if (ngood != 0) { 
                  listb[(addb+1):(addb+ngood),1] <- bonds[w,1]
                  bonds[w,] <- (-1)*(nclust+1)
                  addb <- addb+ngood;
                }
                donea <- donea + 1
              }
              if (doneb != addb){ 
                w <- which(bonds[,1] == listb[doneb+1])
                ngood <- length(w);
                if (ngood != 0) {
                  lista[(adda+1):(adda+ngood),1] <- bonds[w,2]
                  bonds[w,] <- (-1)*(nclust+1)
                  adda <- adda+ngood;
                }
                doneb <- doneb + 1
              }
              if ((donea != adda) || (doneb != addb)){ 
                true <- FALSE
              } else {  
                true = TRUE
              }
            }
            
            pp <- sort(listb[1:doneb])
            pqx <- order(listb[1:doneb])
            #%unx =  unq(listb(1:doneb),pqx);
            #%implanting unq directly
            arr <- listb[1:doneb]
            q <- arr[pqx]
            indices <- which(q != circshift(q,-1))
            count <- length(indices)
            if (count > 0){
              unx <- pqx[indices]
            } else {
              unx <- length(q) -1
            }
            
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            xsz <- length(unx)
            
            pp <- sort(lista[1:donea])
            pqy <- order(lista[1:donea])
            
            #%uny =  unq(lista(1:donea),pqy);
            #%implanting unq directly
            arr <- lista[1:donea]
            q <- arr[pqy]
            indices <- which(q != circshift(q,-1))
            count <- length(indices)
            if (count > 0){
              uny <- pqy[indices]
            } else {
              uny <- length(q) -1
            }
            
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            ysz <- length(uny)
            if ((xsz*ysz) > maxsz){
              maxsz <- xsz*ysz
              mxsz <- xsz
              mysz <- ysz 
            }
            
            thru <- thru - xsz
            nclust <- nclust + 1
          }
          bmap <- bonds[,2]
        }
        
        #% THE SUBNETWORK CODE ENDS
        #% put verbose in for Jaci
        
        ## Adjusting nclust
        all_clusts <- unique(abs(bmap))
        nclust <- length(all_clusts)
        
        #%   THE PERMUTATION CODE BEGINS
        for (nc in 1:nclust){
          
          #message(paste0("nc=", nc))
          w <- which(bmap == (-1)*(nc))
          
          nbonds <- length(w)
          bonds <- mbonds[w,]
          lensq <- bondlen[w]
          
          pq <- sort(bonds[,1])
          st <- order(bonds[,1])
          #%un = unq(bonds(:,1),st);
          #%implanting unq directly     
          arr <- bonds[,1]
          q <- arr[st]
          indices <- which(q != circshift(q,-1))
          count <- length(indices)
          if (count > 0) {
            un <- st[indices]
          } else {
            un <- length(q) - 1
          }
          
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          uold <- bonds[un,1]
          nold <- length(uold)
          
          #%un = unq(bonds(:,2));
          #%implanting unq directly  
          indices <- which(bonds[, 2] != circshift(bonds[, 2], -1))
          count <- length(indices)
          if (count > 0){
            un <- indices
          } else {  
            un <- length(bonds[,2]) -1
          }
          
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          
          unew <- bonds[un,2]
          nnew <- length(unew)
          
          if (nnew > 5){
            rnsteps <- 1
            for (ii in 1:nnew){
              rnsteps <- rnsteps * length(which(bonds[,2] == unew[ii]))
              if (rnsteps >= 50000 && rnsteps < 200000){
                
                warn_log[[length(warn_log) + 1]]  <- i
                #message('Warning: difficult combinatorics encountered')
              } else if (rnsteps >= 200000 && !force_exec){
                
                #message(paste0("i=", i, "... nc=", nc))
                #message('Excessive Combinitorics LOOK WHAT YOU HAVE DONE TO ME!!!')
                #close(uniquevar);
                
                warn_message(warn_log = warn_log, quiet = quiet)
                message(paste0("Excessive Combinitorics encountered while processing slide ", i, 
                               ". Quitting now... Try using a smaller maxdisp."))
                
                return(NULL)
                
              } else if (rnsteps < 5000000 && force_exec) {
                
                warn_log[[length(warn_log) + 1]]  <- i
                
              } else if (rnsteps >= 200000) {
                
                warn_message(warn_log = warn_log, quiet = quiet)
                message(paste0("Excessive Combinitorics encountered while processing slide ", i, 
                               ". Quitting now... Try using a smaller maxdisp."))
                
                return(NULL)
              }
            }
          }
          
          st <- rep(0, times = nnew)
          fi <- rep(0, times = nnew)
          
          h <- rep(0, times = nbonds)
          ok <- rep(1, times = nold)
          nlost <- (nnew - nold) > 0
          
          
          for (ii in 1:nold) { 
            h[which(bonds[,1] == uold[ii])] <- ii
          }
          st[1] <- 1
          fi[nnew] <- nbonds; ##----------------------% check this later
          if (nnew > 1){ 
            sb <- bonds[, 2]
            sbr <- circshift(sb,1)
            sbl <- circshift(sb,-1)
            st[2:length(st)] <- which(sb[2:length(sb)] != sbr[2:length(sbr)]) + 1
            fi[1:(nnew-1)] <- which(sb[1:(nbonds-1)] != sbl[1:(nbonds-1)])
          }
          #%                if i-1 == 13
          #%                    hi
          #%                end
          
          
          checkflag <- 0
          while (checkflag != 2){
            
            pt <- st - 1
            lost <- matrix(0, nrow = nnew, ncol = 1)
            who <- 0
            losttot <- 0
            mndisq <- nnew*maxdisq
            
            
            while (who != (-1)){
              
              if (pt[(who+1)] != fi[(who+1)]){
                
                
                w <- which(ok[h[(pt[(who+1)]+1):(fi[(who+1)])]]!=0) ###---------------% check this -1
                ngood <- length(w)
                if (ngood > 0){
                  if (pt[(who+1)] != st[(who+1)]-1) {
                    ok[h[pt[(who+1)]]] <- 1
                  }
                  pt[(who+1)] <- pt[(who+1)] + w[1]
                  ok[h[pt[(who+1)]]] <- 0
                  if (who == (nnew - 1)){
                    ww <- which(lost == 0)
                    dsq <- sum(lensq[pt[ww]]) + (losttot * maxdisq)
                    
                    if (dsq < mndisq){ 
                      minbonds <- pt[ww]
                      mndisq <- dsq
                    }
                  } else {
                    who <- who+1
                  }
                } else {
                  if (!lost[(who+1)] & (losttot != nlost)){
                    lost[(who+1)] <- 1
                    losttot <- losttot + 1
                    if (pt[(who+1)] != st[(who+1)] - 1){
                      ok[h[pt[(who+1)]]] <- 1
                    }
                    if (who == (nnew-1)){
                      ww <- which(lost == 0)
                      dsq <- sum(lensq[pt[ww]]) + (losttot * maxdisq)
                      if (dsq < mndisq){
                        minbonds <- pt[ww]
                        mndisq <- dsq
                      }
                    } else {
                      who <- who + 1
                    }
                    
                  } else {
                    if (pt[(who+1)] != (st[(who+1)] - 1)){ 
                      ok[h[pt[(who+1)]]] <- 1
                    }
                    pt[(who+1)] <- st[(who+1)] - 1
                    if (lost[(who+1)]){ 
                      lost[(who+1)] <- 0
                      losttot <- losttot -1
                    }
                    who <- who - 1
                  }
                }
                
              } else {  
                if (!lost[(who+1)] && (losttot != nlost)){
                  lost[(who+1)] <- 1
                  losttot <- losttot + 1
                  if (pt[(who+1)] != st[(who+1)]-1){
                    ok[h[pt[(who+1)]]] <- 1
                  }
                  if (who == (nnew - 1)) {
                    ww <- which(lost == 0)
                    dsq <- sum(lensq[pt[ww]]) + (losttot * maxdisq)
                    
                    if (dsq < mndisq){
                      minbonds <- pt[ww]
                      mndisq <- dsq
                    }
                  } else {
                    who <- who + 1
                  }
                } else {
                  if (pt[(who+1)] != st[(who+1)] - 1){
                    ok[h[pt[(who+1)]]] <- 1
                  }
                  pt[(who+1)] <- st[(who+1)] - 1 
                  if (lost[(who+1)]){ 
                    lost[(who+1)] <- 0
                    losttot <- losttot - 1
                  }
                  who <- who -1
                }
              }
            }
            
            checkflag <- checkflag + 1
            if (checkflag == 1){
              plost <- min(c(matfix(mndisq/maxdisq) , (nnew -1)))
              if (plost > nlost){ 
                nlost <- plost 
              } else {
                checkflag <- 2
              }
            }
          }  
          #%   update resx using the minimum bond configuration               
          
          resx[ispan, labely[bonds[minbonds, 2]]] <- eyes[labelx[(bonds[minbonds,1] + 1)]]
          found[labelx[(bonds[minbonds,1] + 1)], 1] <- 1
          
        }
        
        #%   THE PERMUTATION CODE ENDS
      }
      
      w <- which(resx[ispan,] >= 0)
      nww <- length(w)
      
      if (nww > 0){ 
        pos[w,] <- xyzs[resx[ispan,w], (1:dim)]
        if (goodenough > 0){ 
          nvalid[w] <- nvalid[w] + 1
        }
      }  #-----------------------  %go back and add goodenough keyword thing   
      newguys <- which(found == 0)
      nnew <- length(newguys)
      
      if (nnew > 0) {             ##% & another keyword to workout inipos
        newarr <- matrix(-1, nrow = zspan, ncol = nnew)
        
        # cbind?
        resx <- cbind(resx, newarr)
        
        resx[ispan, ((n+1):ncol(resx))] <- eyes[newguys]
        pos <- rbind(pos, xyzs[eyes[newguys],(1:dim)])
        nmem <- matrix(0, nrow = nnew, ncol = 1)
        mem <- c(mem, nmem)
        nun <- 1:nnew
        uniqid <- c(uniqid, ((nun) + maxid))
        maxid <- maxid + nnew
        if (goodenough > 0){ 
          dumphash <- c(dumphash, t(matrix(0, nrow = 1, ncol = nnew)))
          nvalid <- c(nvalid, t(matrix(1, nrow = 1, ncol = nnew)))
        }
        
        #% put in goodenough 
        n <- n + nnew
        
      }
      
    } else {
      #' Warning- No positions found for t='
      message("@@", appendLF = FALSE)
    }
    
    w <- which(resx[ispan,] != (-1))
    nok <- length(w)
    if (nok != 0){
      mem[w] <- 0
    }
    
    #---------------------------------------------------
    mem <- mem + (0 + (cbind(resx[ispan,]) == -1))
    wlost <- which(mem == memory_b+1)
    nlost <- length(wlost)
    
    if (nlost > 0){
      pos[wlost, ] <- (-maxdisp)
      if (goodenough > 0){
        wdump <- which(nvalid[wlost] < goodenough)
        ndump <- length(wdump);
        if (ndump > 0){
          dumphash[wlost[wdump]] <- 1
        }
      }
      #% put in goodenough keyword stuff if 
    }
    
    if ((ispan == zspan) | (i == z)){
      nold <- length(bigresx[1,])
      nnew <- n - nold;
      if (nnew > 0){
        
        newarr <- matrix(-1, nrow = z, ncol = nnew)
        
        ## bigresx <- c(bigresx, newarr)
        bigresx <- cbind(bigresx, newarr)
      }
      if (goodenough > 0){  
        if ((sum(dumphash)) > 0){
          wkeep <- which(dumphash == 0)
          nkeep <- length(wkeep)
          resx <- resx[ ,wkeep]
          bigresx <- bigresx[, wkeep]
          pos <- pos[wkeep, ]
          mem <- mem[wkeep]
          uniqid <- uniqid[wkeep]
          nvalid <- nvalid[wkeep]
          n <- nkeep
          dumphash <- matrix(0, nrow = nkeep, ncol = 1)
        }
      }
      
      #% again goodenough keyword
      if (!quiet) {
        
        message(paste0(i, ' of ' , z, ' done. Tracking ', ntrk, ' particles. ', n, ' tracks total.'))
        
      }
      
      if (!is.matrix(bigresx) || nrow(resx) > nrow(bigresx)) {
        bigresx <- rbind(bigresx)
        bigresx <- rbind(bigresx, 
                         matrix(-1, nrow = (nrow(resx) - nrow(bigresx)), 
                                ncol = ncol(bigresx)))
      }
      
      bigresx[(i-(ispan)+1):i,]  <- resx[1:ispan,]
      resx <- matrix((-1), nrow = zspan, ncol = n)
      
      wpull <- which(pos[ ,1] == (-1 * maxdisp))
      npull <- length(wpull)
      
      if (npull > 0){
        lillist <- list()
        for (ipull in 1:npull){
          wpull2 <- which(bigresx[, wpull[ipull]] != (-1))
          npull2 <- length(wpull2)
          
          
          thing = cbind(bigresx[wpull2,wpull[ipull]],
                        rep(x = uniqid[wpull[ipull]], times = npull2))
          
          #thing <- c(bigresx[wpull2, wpull[ipull]], ),zeros(npull2,1)+uniqid(wpull(ipull))];
          #lillist = [lillist;thing];
          lillist[[length(lillist) + 1]] <- thing
          
        }
        olist[[length(olist) + 1]] <- do.call(rbind, lillist) 
        
      }
      
      
      
      wkeep <- which(pos[, 1] >= 0)
      nkeep <- length(wkeep)
      if (nkeep == 0) { 
        message ('Were going to crash now, no particles....')
      }
      resx <- resx[,wkeep]
      bigresx <- bigresx[, wkeep]
      pos <- pos[wkeep, ]
      mem <- mem[wkeep]
      uniqid <- uniqid[wkeep]
      n <- nkeep;
      dumphash <- matrix(0, nrow = nkeep, ncol =1)
      if (goodenough > 0){
        nvalid <- nvalid[wkeep]
      }
    }
    #waitbar(i / z)  
  }            
  
  if (goodenough > 0){ 
    nvalid <- apply(bigresx >= 0 , 2, sum)
    wkeep <- which(nvalid >= goodenough)
    nkeep <- length(wkeep)
    if (nkeep == 0){
      for (i in 1:10){
        message('You are not going any further, check your params and data')
      }
      message('the code broke at line 1995') 
      return()
    }
    if (nkeep < n){
      bigresx <- bigresx[, wkeep]
      n <- nkeep
      uniqid <- uniqid[wkeep]
      pos <- pos[wkeep, ]
    }
  }
  
  wpull <- which(pos[, 1] != ((-2) * maxdisp))
  npull <- length(wpull);
  if (npull > 0) {
    lillist <- list()
    for (ipull in 1:npull){
      wpull2 <- which(bigresx[, wpull[ipull]] != (-1))
      npull2 <- length(wpull2)   
      thing <- cbind(bigresx[wpull2, wpull[ipull]], 
                     rep(uniqid[wpull[ipull]], times = npull2)) 
      lillist[[length(lillist) + 1]] <- thing
    }
    
    olist[[length(olist) + 1]] <- do.call(rbind, lillist)
  }
  
  olist <- do.call(rbind, olist)
  #%bigresx = 0;
  #%resx = 0;
  
  nolist <- nrow(olist)
  res <- matrix(0, nrow = nolist, ncol = (dd+1))
  for (j in 1:dd){
    res[, j] <- xyzs[olist[, 1], j]
  }
  res[, dd+1] <- olist[,2]
  
  #% this is uberize included for simplicity of a single monolithic code
  
  ndat <- ncol(res)
  newtracks <- res
  
  
  #%u=unq(newtracks(:,ndat));
  
  #% inserting unq
  indices <- which(newtracks[, ndat] != circshift(newtracks[, ndat], -1))
  count <- length(indices)
  if (count > 0){
    u <- indices
  } else {  
    u <- nrow(newtracks) - 1
  }
  
  
  ntracks <- length(u)
  u <- c(0, u)
  for (i in 2: (ntracks + 1)){
    newtracks[(u[(i-1)]+1):u[i], ndat] = (i - 1)
  }
  
  #% end of uberize code
  warn_message(warn_log = warn_log, quiet = quiet)
  return(newtracks)
}


#' Compute Cell Migration Statistics
#'
#' Calculate the statistics from X/Y positional data obtained from cell tracks
#'
#' @param tracks data.frame with cell tracks information
#' @param interval_time integer, time interval between two successive frames were taken
#' @param pixel_micron integer, image resolution, i.e. number of pixels per micron
#'
#' @return list of stats calculated for the cell tracks. Info include variables of speed, 
#' distance, euclidean displacement, persistence, angular displacement, 
#' yFMI, xFMI, y-displacement, x-displacement and frames
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @keywords internal
migrationStats <- function(tracks, interval_time, pixel_micron) {
  
  #
  speed <- list()
  distance <- list()
  frames <- list()
  euclid <- list()
  persistence <- list()
  initial <- list()
  final <- list()
  deltaX = list()
  deltaY = list()
  
  # Keep track
  kept <- list()
  
  cell_number <- sort(unique(tracks[,4]))
  
  for (i in cell_number){
    data <- tracks[tracks[,4] == i, ]
    
    if (!"matrix" %in% class(data) ) {
      next()
    } else (
      kept[[length(kept) + 1]] <- i
    )
    
    # obtain X-axis and Y-axis positional data for the specified track
    X <- data[,1] 
    Y <- data[,2]
    
    x1 <- X[1]
    xEnd <- X[length(X)]
    y1 <- Y[1]
    yEnd <- Y[length(Y)]
    
    initial[[length(initial) + 1]] <- c(x=x1, y=y1)
    final[[length(final) + 1]] <- c(x=xEnd, y=yEnd)
    delX <- as.numeric(xEnd - x1)
    delY <- as.numeric(yEnd - y1)
    deltaX[[length(deltaX) + 1]] <- delX
    deltaY[[length(deltaY) + 1]] <- delY
    # calculate euclidean distance (vector displacement of the cell)
    E <- as.numeric(sqrt((delX)^2 + (delY)^2))
    euclid[[length(euclid) + 1]] <- E
    
    # add subsequent displacements of the cell
    cumulative_displacements <- as.numeric(cumsum(sqrt(diff(X)^2 + diff(Y)^2)) )
    
    # sum of the displacements between each cell centroid for the given
    # track
    
    distance[[length(distance) + 1]] <- max(cumulative_displacements) 
    
    # calculate cell persistence
    persistence[[length(persistence) + 1]] = E/max(cumulative_displacements) 
    
    # total number of frames that cell centroid was tracked ( can be
    # greater than number of frames where centroid was identified given the
    # param.mem parameter
    # total number of time intervals through which the cell has been tracked
    totalframes = data[nrow(data), 3] - data[1, 3] 
    
    # sum of all individual displacemnts divided by the time that cell
    # centroid was tracked
    ds_dt <- max(cumulative_displacements)/(totalframes*interval_time)
    speed[[length(speed) + 1]] <- ds_dt
    
    frames[[length(frames) + 1]] <- totalframes
  }
  
  # Expand resulting lists
  speed <- do.call(c, speed)
  distance <- do.call(c, distance)
  frames <- do.call(c, frames)
  euclid <- do.call(c, euclid)
  persistence <- do.call(c, persistence)
  initial <- do.call(rbind, initial)
  final <- do.call(rbind, final)
  kept <- do.call(c, kept)
  deltaX <- do.call(c, deltaX)
  deltaY <- do.call(c, deltaY)
  
  
  # calculate angular displacement of cells trajectory
  arccos <- deltaX/euclid
  theta <- acos(arccos)
  
  for (j in 1:length(arccos)){
    if  (arccos[j] < 0 & deltaY[j] > 0) {
      theta[j] <- 2*pi - theta[j]      
    }
    
    if (arccos[j] > 0 & deltaY[j] > 0) {
      theta[j] <- 2*pi - theta[j]      
    }
  }
  
  #theta = theta.*((2*pi)/360);
  deltaY <-  deltaY * (-1)
  yfmi <- deltaY / distance
  xfmi <- deltaX / distance
  deltaY <- deltaY * pixel_micron
  deltaX <- deltaX * pixel_micron
  speed <- speed * pixel_micron;
  distance <- distance * pixel_micron
  euclid <- euclid * pixel_micron
  
  OUT <- list(
    speed = speed,
    distance = distance,
    frames = frames,
    euclid = euclid,
    persistence = persistence,
    initial = initial,
    final = final,
    yfmi = yfmi,
    xfmi = xfmi,
    deltaX = deltaX, 
    deltaY = deltaY, 
    kept = kept,
    theta = theta
  )
  return(OUT)
}






#' Compute Tracks Stats
#'
#' Wrapper for the migrationStats() function. It computes statistics for a 
#' trackedCells object where cells have already been tracked.
#'
#' @param tc_obj a trackedCells object
#' @param time_between_frames integer, time interval between two successive frames were taken
#' @param resolution_pixel_per_micron integer, image resolution, i.e. number of pixels per micron
#' 
#' @return a trackedCells object, including cell track statistics 
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @export
compute_tracks_stats <- function(tc_obj, time_between_frames, resolution_pixel_per_micron)
{
  
  if(tc_obj@ops$track == 0)
    stop("You need to run cell_tracker() before computing stats")
  
  # RETRIEVE
  my_tracks <- tc_obj@tracks
  
  # DO
  handles <- migrationStats(tracks = my_tracks, interval_time = time_between_frames, 
                            pixel_micron = resolution_pixel_per_micron);
  
  
  sz <- length(handles$speed)
  handles$cell_number <- 1:sz
  
  cell_stats <- data.frame(Cell_Number = handles$cell_number, 
                           Speed = handles$speed, 
                           Distance = handles$distance, 
                           Displacement = handles$euclid, 
                           Persistence = handles$persistence, 
                           Degrees = handles$theta, 
                           YFMI = handles$yfmi, 
                           XFMI = handles$xfmi, 
                           Y_displacement = handles$deltaY, 
                           X_displacement = handles$deltaX, 
                           Frames = handles$frames, 
                           stringsAsFactors = FALSE)
  
  # to organize population stats
  my_colz <- c("Speed", "Distance", "Displacement", "Persistence", "YFMI", 
               "XFMI", "Y_displacement", "X_displacement")
  
  my_rows <- lapply(my_colz, function(cl) {
    tmp <- cell_stats[, cl]
    data.frame(mean = mean(tmp, na.rm = TRUE),
               SD = sd(tmp, na.rm = TRUE), 
               median = median(tmp, na.rm = TRUE), 
               min = min(tmp, na.rm = TRUE), 
               max = max(tmp, na.rm = TRUE))
  })
  my_rows <- do.call(rbind, my_rows)
  rownames(my_rows) <- my_colz
  
  # compute sum of cos and sin of angles
  r <- sum(exp(1i*handles$theta))
  
  # obtain mean angle 
  meanTheta <- Arg(r)
  degrees <- meanTheta/pi*180
  my_rows <- rbind(my_rows, 
                   Angle = data.frame(mean = sz, SD = degrees, median = NA, min = NA, max = NA))
  rownames(my_rows)[c(2,3)] <- c("Total_displacement", "Euclidean_displacement")
  pop_stats <- my_rows
  
  # Attach, return
  tc_obj@ops$stats <- 1
  tc_obj@stats <- list(population = pop_stats, 
                       cells = cell_stats)
  
  return(tc_obj)
}





# --------------------------------------------------
# Module #5 - User interface
# --------------------------------------------------


#' Optimize Detection Params
#'
#' Optimize Detection Parameters for running a cell tracking job
#'
#' @details The lnoise param is used to guide a lowpass blurring operation, while the lobject param is used
#' to guide a highpass background subtraction. The threshold param is used for a background correction following
#' the initial image convolution
#' \itemize{
#' 
#'   \item \strong{lnoise}: Characteristic lengthscale of noise in pixels. 
#'   Additive noise averaged over this length should vanish. May assume any positive floating value.
#'   May be also set to 0, in which case only the highpass "background subtraction" operation is performed.
#' 
#'   \item \strong{lobject} Integer length in pixels somewhat larger than a typical object. 
#'   Can also be set to 0, in which case only the lowpass "blurring" operation defined by lnoise is done 
#'   without the background subtraction defined by lobject
#'   
#'   \item \strong{threshold} Numeric. By default, after the convolution, any negative pixels are reset 
#'   to 0.  Threshold changes the threshhold for setting pixels to 0.  Positive values may be useful 
#'   for removing stray noise or small particles. 
#'   
#' }
#' 
#'
#'
#' @param tc_obj a trackedCells object
#' @param lnoise_range numeric vector of lnoise values to be used in the optimization step. Can be NULL
#' @param diameter_range numeric vector of diameter values to be used in the optimization step. Can be NULL
#' @param threshold_range numeric vector of threshold values to be used in the optimization step. Can be NULL
#' @param target_cell_num integer, the expected (optimal) number of cells to be detected in each frame
#' @param threads integer, number of cores to use for parallelization
#' @param quantile.val numeric, argument passed to estimate_diameter_range(). If NULL, it is defaulted to 0.99 
#' @param px.margin numeric, argument passed to estimate_diameter_range(). If NULL, it ia defaulted to 2

#'
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @importFrom stats quantile
#' @importFrom utils head
#' @importFrom graphics par
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#' @importFrom doParallel registerDoParallel  
#' @import foreach
#'
#' @export
optimize_params <- function(tc_obj, lnoise_range = NULL, 
                            diameter_range = NULL, threshold_range = NULL, 
                            target_cell_num = NULL, threads = 1,
                            quantile.val = NULL, px.margin= NULL)
  
{
  # do
  stack_img <- tc_obj@images
  
  # Nested f(x)
  all_combos <- function(...){
    xx <- list(...)
    zz <- names(xx)
    
    # Init
    out <- data.frame(xx[[1]], stringsAsFactors = FALSE)
    colnames(out) <- zz[1]
    
    # Keep attaching
    for (j in 2:length(xx)) {
      TMP <- xx[[j]]
      nuOUT <- list()
      for (z in TMP){
        tmpout <- cbind(out, data.frame(TMP = z, stringsAsFactors = FALSE))
        nuOUT[[length(nuOUT) + 1]] <- tmpout
      }
      out <- do.call(rbind, nuOUT)
      colnames(out)[j] <- zz[j]
    }
    return(out)  
  }
  
  ## ----- debugging -----
  #bpass = cellTracker:::bpass
  #pkfnd = cellTracker:::pkfnd
  #visualize_img = cellTracker:::visualize_img
  #cntrd = cellTracker:::cntrd
  #next_odd = cellTracker:::next_odd
  #visualize_cntr = cellTracker:::visualize_cntr
  #track = cellTracker:::track
  ## ----- endo of debugging -----
  
  # select mid signal image
  imgSums <- sapply(stack_img$images, sum, na.rm = TRUE)
  med.i <- ifelse(length(imgSums) %% 2 == 0, length(imgSums) / 2, (0.5 * length(imgSums) + 0.5))
  r.i <- order(imgSums)[med.i]
  
  tmp_img <- stack_img$images[[r.i]]
  
  # Define param ranges
  if (is.null(px.margin)) {
    px.margin <- 2
  }
  if (is.null(quantile.val)) {
    quantile.val <- 0.99
  }
  
  estRDI <- tryCatch({
    estimate_diameter_range(x = tmp_img, px.margin = px.margin,
                            quantile.val = quantile.val, plot = FALSE)}, 
    error = function(e) NULL)
  
  # diam range
  if(is.null(diameter_range) && !is.null(estRDI)) {
    diameter_range <- c(floor(estRDI$q50.diam - 1), ceiling(estRDI$q95.diam + 5))
    diameter_range[diameter_range < 2] <- 2
    diameter_range <- unique(as.integer(diameter_range))
    diameter_range <- unique(as.integer(
      seq(min(diameter_range), max(diameter_range), length.out = 3))) 
    
  } else if (is.null(diameter_range)) {
    diameter_range <- c(10, 50, 100)
  }
  
  # num cell
  if(is.null(target_cell_num) && !is.null(estRDI)) {
    target_cell_num <- estRDI$estim.cell.num
  } else if (is.null(target_cell_num)) {
    target_cell_num <- 100
  }
  
  # Define param ranges
  if(is.null(lnoise_range))
    lnoise_range <- c(2, 8, 16)
  

  if(is.null(threshold_range)) {
    threshold_range <- seq(max(0, (min(tmp_img[tmp_img > min(tmp_img, na.rm = TRUE)], na.rm = TRUE) - 1)), 
                           (1 + quantile(tmp_img[tmp_img > min(tmp_img, na.rm = TRUE)], probs = 0.75)),
                           length.out = 4)
    threshold_range <- unique(as.integer(threshold_range))
  }
  
  # Al params
  all_params <- all_combos(image.i = med.i,
                           lnoise = lnoise_range, 
                           diameter = diameter_range,
                           threshold = threshold_range)
  
  all_params <- all_params[all_params$diameter > (4 + all_params$lnoise), ]
  rownames(all_params) <- NULL
  
  if(nrow(all_params) < 1) {
    message("There is a problem with the param ranges that were submitted")
    message("Please, try again with different param ranges")
    return(tc_obj)
  }
  
  # Verbose
  message(paste0("Testing ", nrow(all_params), " combination(s) of params."), appendLF = TRUE)
  message("This may take some time.", appendLF = TRUE)
  
  message("Processing ", appendLF = FALSE)
  
  ##
  ## Parallelize please
  j <- NULL
  
  # how many cores can we use?
  num_parallelCores <- threads
  debugging <- TRUE
  
  max.cores <- parallel::detectCores()
  max.cores <- max.cores - 1
  max.cores <- ifelse(max.cores < 1, 1, max.cores)
  my.test <- 1 <= num_parallelCores & num_parallelCores <= max.cores
  use.cores <- ifelse(my.test, num_parallelCores, max.cores)
  
  # cores = 1, do not parallelize
  if (use.cores == 1) {
    
    # Initialize collector (list)
    all_results <- list()
    
    for (i in 1:nrow(all_params)){
      
      # Verbose
      message(".", appendLF = FALSE)
      
      #visualize_img(tmp_img)
      b <- bpass(image_array = tmp_img, 
                 lnoise = all_params$lnoise[i], 
                 lobject = all_params$diameter[i], 
                 threshold = all_params$threshold[i])
      tmpOUT <- list(img = b)
      
      tryCatch({
        pk <- suppressMessages(
          pkfnd(im = b, 
                th = all_params$threshold[i], 
                sz = next_odd(all_params$diameter[i])))
        
        cnt <- suppressMessages(
          cntrd(im = b, mx = pk, 
                sz = next_odd(all_params$diameter[i])))
        
        tmpOUT[["count"]] <- nrow(cnt)
        
      }, error = function(e) {
        tmpOUT[["count"]] <- 0
        
      })
      all_results[[i]] <- tmpOUT
    }
    
    # cores > 1, DO parallelize!
  } else {
    
    if (debugging) {
      cl <- suppressMessages(parallel::makeCluster(use.cores, outfile = ""))
    } else {
      cl <- suppressMessages(parallel::makeCluster(use.cores))
    }
    
    suppressMessages(doParallel::registerDoParallel(cl))
    
    # Nothing to export! "tmp_img", "all_params" automatically exported
    #stuffToExp <- c("tmp_img", "all_params")
    stuffToExp <- c()
    suppressMessages(parallel::clusterExport(cl, stuffToExp))
    
    ## %dopar%
    all_results <- 
      tryCatch(foreach::foreach(j = (1:nrow(all_params)),
                                .verbose = TRUE, 
                                .packages = "cellTracker") %dopar% {
                                  
                                  # Verbose
                                  message(".", appendLF = FALSE)
                                  
                                  #visualize_img(tmp_img)
                                  b <- bpass(image_array = tmp_img, 
                                             lnoise = all_params$lnoise[j], 
                                             lobject = all_params$diameter[j], 
                                             threshold = all_params$threshold[j])
                                  tmpOUT <- list(img = b)
                                  
                                  tryCatch({
                                    pk <- suppressMessages(
                                      pkfnd(im = b, 
                                            th = all_params$threshold[j], 
                                            sz = next_odd(all_params$diameter[j])))
                                    
                                    cnt <- suppressMessages(
                                      cntrd(im = b, mx = pk, 
                                            sz = next_odd(all_params$diameter[j])))
                                    
                                    tmpOUT[["count"]] <- nrow(cnt)
                                    
                                  }, error = function(e) {
                                    tmpOUT[["count"]] <- 0
                                    
                                  })
                                  tmpOUT
                                }, error = (function(e) {
                                  print(e)
                                  try(parallel::stopCluster(cl), silent = TRUE)
                                  return(NULL)
                                }))
    message("Done!", appendLF = TRUE)
    try({suppressWarnings(parallel::stopCluster(cl))}, silent = TRUE)
  }
  
  # Attach counts
  all_params$counts <- do.call(c, lapply(all_results, function(x) {
    tmp <- x$count
    ifelse(is.null(tmp), 0, tmp)}))
  all_params$i <- 1:nrow(all_params)
  
  # Return
  all_params$diff100 <- abs(target_cell_num - all_params$counts)
  ord_params <- all_params[order(all_params$diff100), ]
  ret.i <- head(ord_params$i, n = 9)
  best_params <- list()
  
  curPAR <- par(no.readonly = TRUE)
  par(mfrow = c(3, 3))
  top.i <- 1
  for (ri in ret.i) {
    
    if (top.i == 1) {
      best_params[["lnoise"]] <- ord_params$lnoise[ord_params$i == ri]
      best_params[["diameter"]] <- ord_params$diameter[ord_params$i == ri]
      best_params[["threshold"]] <- ord_params$threshold[ord_params$i == ri]
    }
    
    myLAB <- paste0("Pick #", top.i, "; Cell_count=", ord_params$counts[ord_params$i == ri], "\n")
    myLAB <- paste0(myLAB, "lnoise=", ord_params$lnoise[ord_params$i == ri], "; ",
                    "diameter=", ord_params$diameter[ord_params$i == ri], "; ",
                    "threshold=", ord_params$threshold[ord_params$i == ri])
    
    visualize_img(img_mtx = all_results[[ri]]$img, 
                  main = myLAB)
    
    top.i <- top.i + 1
  }
  par(curPAR)
  
  # Extract_all_img 
  allIMG <- lapply(all_results, function(x) {x$img})
  
  #return(list(auto_params = best_params, 
  #            results = all_params, 
  #            images = allIMG))
  
  
  tc_obj@ops$optimized_params <- 1
  tc_obj@optimized <- list(auto_params = best_params, 
                           results = all_params)
  
  return(tc_obj)
}


#' Compute Cell Tracks 
#'
#' Analyze Stacks, detect cells in each frame, and analyze cell tracks over time
#'
#' @details The lnoise param is used to guide a lowpass blurring operation, while the lobject param is used
#' to guide a highpass background subtraction. The threshold param is used for a background correction following
#' the initial image convolution
#' \itemize{
#' 
#'   \item \strong{lnoise}: Characteristic lengthscale of noise in pixels. 
#'   Additive noise averaged over this length should vanish. May assume any positive floating value.
#'   May be also set to 0, in which case only the highpass "background subtraction" operation is performed.
#' 
#'   \item \strong{lobject} Integer length in pixels somewhat larger than a typical object. 
#'   Can also be set to 0, in which case only the lowpass "blurring" operation defined by lnoise is done 
#'   without the background subtraction defined by lobject
#'   
#'   \item \strong{threshold} Numeric. By default, after the convolution, any negative pixels are reset 
#'   to 0.  Threshold changes the threshhold for setting pixels to 0.  Positive values may be useful 
#'   for removing stray noise or small particles. 
#'   
#' }
#' 
#'
#' @param tc_obj a trackedCells object
#' @param lnoise numeric, lnoise parameter; can be NULL if optimize_params() has already been run
#' @param diameter numeric, diameter parameter; can be NULL if optimize_params() has already been run
#' @param threshold numeric, threshold parameter; can be NULL if optimize_params() has already been run
#' @param maxDisp numeric,  maximum displacement of a cell per time interval. 
#' When many cells are detected in each frame, small maxDisp values should be used. 
#' @param memory_b numeric, memory_b parameter as used in the original track.m function.
#' In the current R implementation, only the value memory_b=0 is accepted
#' @param goodenough numeric, goodenough parameter as used in the original track.m function.
#' In the current R implementation, only the value goodenough=0 is accepted
#' @param threads integer, number of cores to use for parallelization
#' @param show_plots logical, shall cells detected in each frame of the image stack be visualized
#' @param verbose logical, shall info about the progress of the cell tracking job be printed
#' 
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#' @importFrom doParallel registerDoParallel  
#' @import foreach
#'
#' @export
cell_tracker <- function(tc_obj, 
                         lnoise = NULL, diameter = NULL, 
                         threshold = NULL, maxDisp = 25, 
                         memory_b = 0, goodenough = 0, 
                         threads = 1,
                         show_plots = TRUE, verbose = TRUE) 
{
  # get stuff
  stack_img <- tc_obj@images
  optimal_params <- tc_obj@optimized
  
  if (length(optimal_params) > 0) {
    my.lnoise <- optimal_params$auto_params$lnoise
    my.diameter <- optimal_params$auto_params$diameter
    my.threshold <- optimal_params$auto_params$threshold
  } else {
    my.lnoise <- NA
    my.diameter <- NA
    my.threshold <- NA
  }
  
  custom_params_flag <- 0
  
  if(!is.null(lnoise) && is.numeric(lnoise)){
    my.lnoise <- lnoise[1]
    custom_params_flag <- 1
  }
  
  if(!is.null(diameter) && is.numeric(diameter)){
    my.diameter <- diameter[1]
    custom_params_flag <- 1
  }
  
  if(!is.null(threshold) && is.numeric(threshold)){
    my.threshold <- threshold[1]
    custom_params_flag <- 1
  }
  
  # Max Disp
  if(!is.null(maxDisp) && is.numeric(maxDisp)) {
    maxDisp <- maxDisp[1]
  } else {
    maxDisp <- 20
  }
  
  # Other params
  quiet <- verbose 
  force_exec <- FALSE
  j <- NULL
  
  # At this moment, let's play safe!
  # Impose the following
  if (memory_b != 0)
    message("Currently, only memory_b=0 is supported... Resetting.")
  memory_b <- 0 
  
  if (goodenough != 0)
    message("Currently, only goodenough=0 is supported... Resetting.")
  goodenough <- 0
  
  # In the end, my params are:
  lnoise <- my.lnoise
  diameter <- my.diameter
  threshold <- my.threshold
  
  track_params <- list(lnoise = lnoise, 
                       diameter = diameter, 
                       threshold = threshold, 
                       maxDisp = maxDisp,
                       goodenough = goodenough, 
                       memory_b = memory_b,
                       force_exec = force_exec,
                       quiet = quiet, 
                       verbose = verbose,
                       show_plots = show_plots)
  
  # Final check
  if (sum(sapply(track_params, is.na)) > 0) {
    message("Make sure to set all params for the analysis, or run optimize_params()")
    return(tc_obj)
  }
  
  ## ----- debugging -----
  #bpass = cellTracker:::bpass
  #pkfnd = cellTracker:::pkfnd
  #visualize_img = cellTracker:::visualize_img
  #cntrd = cellTracker:::cntrd
  #next_odd = cellTracker:::next_odd
  #visualize_cntr = cellTracker:::visualize_cntr
  #track = cellTracker:::track
  ## ----- endo of debugging -----
  
  # Load stack
  stack <- stack_img
  
  InfoImage <- stack$attributes[[1]]
  mImage <- stack$dim$width_m
  nImage <- stack$dim$height_n
  NumberImages <- stack$dim$NumberImages
  FinalImage <- stack$images
  
  ## ----------- Evaluate centroids ---------------------
  
  # locate centroids, via centroid_array
  if (verbose)
    message("Computing centroid positions", appendLF = FALSE)
  
  ##
  ## Parallelize please
  
  # how many cores can we use?
  num_parallelCores <- threads
  debugging <- TRUE
  
  max.cores <- parallel::detectCores()
  max.cores <- max.cores - 1
  max.cores <- ifelse(max.cores < 1, 1, max.cores)
  my.test <- 1 <= num_parallelCores & num_parallelCores <= max.cores
  use.cores <- ifelse(my.test, num_parallelCores, max.cores)
  
  # cores = 1, do not parallelize
  if (use.cores == 1) {
    
    # Init collectors
    all_centroids <- list()
    all_b <- list()
    
    for (i in 1:NumberImages) {
      
      if (verbose)
        message(".", appendLF = FALSE)
      
      # generate an 1xP array with each column containing centroid output for 
      # individual frames
      a <- FinalImage[[i]]
      b <- bpass(image_array = a, lnoise = lnoise, lobject = diameter, threshold = threshold)
      pk <- pkfnd(im = b, th = threshold, sz = next_odd(diameter))
      cnt <- cntrd(im = b, mx = pk, sz = next_odd(diameter))
      
      if (show_plots) {
        visualize_img(img_mtx = b, las = 1, main = paste0("Stack num. ", i))
        visualize_cntr(centroids = cnt, width_px = ncol(b), height_px = nrow(b)) 
      }
      
      # determine that frame s has at least 1 valid centroid
      if(! is.null(cnt) && is.data.frame(cnt) && nrow(cnt) > 0) {
        all_centroids[[length(all_centroids) + 1]] <- cnt
        all_b[[length(all_b) + 1]] <- b
        
      } else {
        message(paste0('No centroids detectd in frame ', i, ' in the current stack'))
        message('Please, check nuclei validation settings for this image stack.')
      }
    }
    if (verbose)
      message("", appendLF = TRUE)
    
  } else {
    
    if (debugging) {
      cl <- suppressMessages(parallel::makeCluster(use.cores, outfile = ""))
    } else {
      cl <- suppressMessages(parallel::makeCluster(use.cores))
    }
    
    suppressMessages(doParallel::registerDoParallel(cl))
    # Nothing to export! ""FinalImage", "all_params" automatically exported
    #stuffToExp <- c("FinalImage", "all_params")
    stuffToExp <- c()
    suppressMessages(parallel::clusterExport(cl, stuffToExp))
    
    ## %dopar%
    all_results <- 
      tryCatch(foreach::foreach(j = (1:NumberImages),
                                .verbose = TRUE, 
                                .packages = "cellTracker") %dopar% {
                                  
                                  # Verbose
                                  message(".", appendLF = FALSE)
                                  
                                  # generate an 1xP array with each column containing centroid output for 
                                  # individual frames
                                  a <- FinalImage[[j]]
                                  b <- bpass(image_array = a, lnoise = lnoise, lobject = diameter, threshold = threshold)
                                  pk <- pkfnd(im = b, th = threshold, sz = next_odd(diameter))
                                  cnt <- cntrd(im = b, mx = pk, sz = next_odd(diameter))
                                  
                                  # determine that frame s has at least 1 valid centroid
                                  if(! is.null(cnt) && is.data.frame(cnt) && nrow(cnt) > 0) {
                                    tmpOUT <- list(cnt = cnt, b = b, j = j)
                                    
                                  } else {
                                    message(paste0('No centroids detectd in frame ', i, ' in the current stack'))
                                    message('Please, check nuclei validation settings for this image stack.')
                                    tmpOUT <- list(cnt = NA, b = b, j = j)
                                    
                                  }
                                  tmpOUT
                                  
                                }, error = (function(e) {
                                  print(e)
                                  try(parallel::stopCluster(cl), silent = TRUE)
                                  return(NULL)
                                }))
    
    message("Done!", appendLF = TRUE)
    try({suppressWarnings(parallel::stopCluster(cl))}, silent = TRUE)
    
    re.idx <- order(do.call(c, lapply(all_results, function(x) {x$j})))
    all_results <- all_results[re.idx]
    all_centroids <- lapply(all_results, function(x) {x$cnt})
    all_b <- lapply(all_results, function(x) {x$b})
    
    # Visualize if needed
    if (show_plots) {
      for (ii in 1:length(all_results)){
        bii <- all_results[[ii]]$b
        cntii <- all_results[[ii]]$cnt
        visualize_img(img_mtx = bii, las = 1, main = paste0("Stack num. ", ii))
        visualize_cntr(centroids = cntii, width_px = ncol(bii), height_px = nrow(bii)) 
        bii<-NULL; cntii<-NULL
      }
    }
  }  
  
  # Position list (reformated centroid data for track.m input)
  OUT_centroids <- all_centroids
  # Remove columns that contain brightness and sqare of radius of gyration
  # this is the equivalent of position.m function
  # Also, create a frame(tau)label for each array of centroid data
  for (ti in 1:length(all_centroids)) { 
    
    # retain only positional data by removing columns 3 and 4
    all_centroids[[ti]] <- all_centroids[[ti]] [, 1:2]
    all_centroids[[ti]]$tau <- ti
  }
  
  # create a matrix that contains centroid data in sequential order by frame(tau)
  pos <- do.call(rbind, all_centroids)
  
  ## generate tracks 
  tracks <- track(xyzs = pos, maxdisp = maxDisp, params = track_params)
  
  # pack and return
  #OUT <- list(images = all_b, 
  #            centroids = OUT_centroids,
  #            positions = pos,
  #            tracks = tracks, 
  #            params = track_params)
  
  tc_obj@proc_images <- list(images = all_b)
  tc_obj@centroids <- OUT_centroids
  tc_obj@positions <- pos 
  tc_obj@tracks <- tracks 
  tc_obj@params <- track_params
  
  tc_obj@ops$track <- 1
  tc_obj@ops$custom_params <- custom_params_flag
  
  return(tc_obj)
}



# --------------------------------------------------
# Module #6 - Getters and setters
# --------------------------------------------------

#' Get Track Data
#'
#' Extract Track Data from a trackedCells object
#'
#' @param tc_obj a trackedCells object
#' @param attach_meta logical, shall metaData be attached to tracks
#'
#' @return a data.frame including cell tracks data
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @export
get_tracks <- function(tc_obj, attach_meta = FALSE) 
{
  tmp <- tc_obj@tracks
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- c("Y", "X", "frame.ID", "cell.ID")
  TMP <- tmp[, c("frame.ID", "X", "Y", "cell.ID")]
  if(attach_meta && nrow(TMP) > 0) {
    
    TMP$tiff_file = tc_obj@metadata$tiff_file
    TMP$experiment = tc_obj@metadata$experiment
    TMP$condition = tc_obj@metadata$condition
    TMP$replicate = tc_obj@metadata$replicate
  } else if (nrow(TMP) < 1) {
    return (NULL)
  }
    
  return(TMP)
}



#' Get Cell population stats
#'
#' Extract cell population statistics from a trackedCells object
#'
#' @param tc_obj a trackedCells object
#'
#' @return data.frame including cell population stats
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @export
get_population_stats <- function(tc_obj) 
{
  if (tc_obj@ops$stats == 1) {
    return(tc_obj@stats$population)  
  } else {
    message("Stats have not been computed yet. Please, run `compute_tracks_stats()`. Thanks.")    
  }
}



#' Get Cell migration stats
#'
#' Extract cell migration statistics from a trackedCells object
#'
#' @param tc_obj a trackedCells object
#'
#' @return data.frame including cell migration stats
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' 
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @export
get_cells_stats <- function(tc_obj) 
{
  if (tc_obj@ops$stats == 1) {
    return(tc_obj@stats$cells)  
  } else {
    message("Stats have not been computed yet. Please, run `compute_tracks_stats()`. Thanks.")    
  }
}


#' Get MetaData
#'
#' Extract MetaData from a trackedCells object
#'
#' @param tc_obj a trackedCells object
#'
#' @return a list including four items: tiff filename, experiment name, condition label, and replicate ID.
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @export
get_cells_meta <- function(tc_obj) 
{
  tmp <- tc_obj@metadata
  return(tmp)
}


#' Get Image Stacks
#'
#' Extract Images Stacks from a trackedCells object
#'
#' @param tc_obj a trackedCells object
#'
#' @return a list including stack images (formatted as numeric matrices)
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @export
get_image_stacks <- function(tc_obj) 
{
  tmp <- tc_obj@images$images
  return(tmp)
}



#' Set MetaData
#'
#' Write/Replace MetaData of a trackedCells object
#'
#' @param tc_obj a trackedCells object
#' @param experiment string, a label to describe the experiment (optional). Can be NULL
#' @param condition string, a label to describe the experimental condition (optional). Can be NULL 
#' @param replicate string, a label to identify the replicate (optional). Can be NULL
#'
#' @return a list including three items: experiment name, condition label, and replicate ID.
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @export
set_cells_meta <- function(tc_obj, experiment = NULL, 
                       condition = NULL, replicate = NULL) 
{
  
  if (is.null(experiment)) {
    experiment <- NA  
  } else {
    experiment <- tryCatch(as.character(experiment[1]), error = function(e) NA)
  }
  
  if (is.null(replicate)) {
    replicate <- NA  
  } else {
    replicate <- tryCatch(as.character(replicate[1]), error = function(e) NA)
  }
  
  if (is.null(condition)) {
    condition <- NA  
  } else {
    condition <- tryCatch(as.character(condition[1]), error = function(e) NA)
  }
  
  FILENM <- tc_obj@metadata$tiff_file
  
  tmp <- list(tiff_file = FILENM, 
              experiment = experiment, 
              condition = condition, 
              replicate = replicate)
  tc_obj@metadata <- tmp
  return(tc_obj)
}


#
# Lates f(x)

#' Aggregate trackedCells Objects
#'
#' Aggregate two or more trackedCells-class objects together. Input objects must carry information of cell
#' tracks (oterwise an error will be raised). All tracks form the different experiments/images are returned in a
#' large data.frame. A new unique ID is assigned to specifically identify each cell track from each image/experiment.
#'
#' @param x a trackedCells-class object where cells have already been tracked
#' @param ... one or more trackedCells-class object(s) where cells have already been tracked
#' @param meta_id_field string, can take one of the following values, c("tiff_file", "experiment", 
#' "condition", "replicate"). Indicates the meta-data column used as unique ID for the image/experiment.
#' Can be abbreviated. Defaults to "tiff_file".
#'
#' @return An aggregate data.frame including all cells that were tracked over two or more images/experiments. 
#' The data.frame includes the following columns: "new.ID", "frame.ID", "X", "Y", "cell.ID", "tiff_name", 
#' "experiment", "condition", "replicate". The "new.ID" uniquely identifies a cell in a given image/experiment.
#'
#' @details each trackedCells-class object passed to this function requires a unique identifier (such as a unique 
#' tiff_file name). Any of the metadata columns can be used as unique ID for an image/experiment. The function
#' will raise an error if non-unique identifiers are found across the input objects.
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @importFrom graphics par image
#'
#' @export
aggregate_trackedCells <- function(x, ..., meta_id_field = c("tiff_file", "experiment", 
                                                             "condition", "replicate")) 
{
  # Inner fxs
  check_trobj <- function(xx) {
    RT <- FALSE
    if (!is.null(xx)) {
      if ("trackedCells" %in% class(xx)) {
        if (nrow(xx@tracks) > 0) {
          RT <- TRUE  
        }
      }  
    }
    return(RT)
  }
  
  compute_mult <- function(xx) {
    zz <- nchar(xx) + 2
    out <- (10 ^ zz)
    return(out)
  }
  
  meta_id_field <- match.arg(arg = meta_id_field, 
                             choices = c("tiff_file", "experiment", 
                                         "condition", "replicate"), 
                             several.ok = FALSE)
  y <- list(...)
  test1 <- FALSE
  if (length(y) > 0) {
    test1 <- sum(do.call(c, lapply(y, check_trobj))) == length(y)
  }
  
  # first check
  stopifnot(check_trobj(x), test1)
  
  big.list <- list(x)
  for(yi in y) {
    big.list[[length(big.list) + 1]] <- yi
  }
  
  # Chek names are different
  all_ids <- do.call(c, lapply(big.list, function(xx) {xx@metadata[[meta_id_field]] }))
  all_ids <- as.character(all_ids)
  unq_ids <- unique(all_ids)
  
  # second check
  stopifnot(length(unq_ids) == length(all_ids))
  
  # Adjust
  my_tracks <- lapply(big.list, get_tracks, attach_meta = TRUE)
  my_tracks <- do.call(rbind, my_tracks)
  my_tracks[,"new.ID"] <- factor(my_tracks[,meta_id_field], levels = unq_ids)
  my_tracks[,"new.ID"] <- as.numeric(my_tracks[,"new.ID"])
  my.mult <- compute_mult(max(my_tracks[, "cell.ID"], na.rm = TRUE))
  my_tracks[,"new.ID"] <- (my.mult * my_tracks[,"new.ID"]) + my_tracks[, "cell.ID"]
  
  keep.colz <- c('new.ID', 'frame.ID', 'X', 'Y', 'cell.ID', 'tiff_file', 'experiment', 'condition', 'replicate')
  out <- my_tracks[, keep.colz]
  rownames(out) <- NULL
  
  return(out)
}


#' Filter an Aggregated Table of Cell Tracks
#'
#' Filter an Aggregated Table (data.frame) of cell tracks (from multiple images/experiments) and 
#' retain cell tracks from images/experiments of interest
#'
#' @param x data.frame, is an aggregated Table of Cell Tracks. Must include the following columns:
#' "new.ID", "frame.ID", "X", "Y", "cell.ID", "tiff_name", "experiment", "condition", "replicate"
#' @param id_list character vector, indicates the IDs (such as tiff_filenames) to be retained in the 
#' output data.frame
#' @param meta_id_field string, can take one of the following values, c("tiff_file", "experiment", 
#' "condition", "replicate"). Indicates the meta-data column used as unique ID for the image/experiment.
#' Can be abbreviated. Defaults to "tiff_file".
#'
#' @return data.frame, a filtered aggregated Table of Cell Tracks  
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#'
#' @examples
#' A <- data.frame(new.ID = 1:10, 
#'                 frame.ID = 10:1, 
#'                 X = sample(1:100, size = 10), 
#'                 Y = sample(1:100, size = 10), 
#'                 cell.ID = c(rep(1, 5), rep(2, 5)), 
#'                 tiff_file= c(rep("ii", 3), rep("jj", 5), rep('kk', 2)))
#' filter_trackedCells(A, id_list = c("ii", "kk"))
#'
#' @export
filter_trackedCells <- function(x, id_list, 
                                meta_id_field = c("tiff_file", "experiment", 
                                                  "condition", "replicate")) {
  
  meta_id_field <- match.arg(arg = meta_id_field, 
                             choices = c("tiff_file", "experiment", 
                                         "condition", "replicate"), 
                             several.ok = FALSE)
  
  REQd <- c("new.ID", "frame.ID", "X", "Y", "cell.ID", meta_id_field)
  CHK1 <- sum(REQd %in% colnames(x)) == length(REQd)
  
  stopifnot(CHK1)
  
  xx <- x[x[, meta_id_field] %in% id_list, ]
  return(xx)
}





#' @title Sample Stack of Fluorescent Cells 
#' @description Sample Stack of Fluorescent Cells to be used for computing cell tracks and stats
#' 
#' @keywords internal
"sample_tracked_cells"


#' @title Track Fluorescent Cells and Compute Migration Stats
#' @description Track Fluorescent Cells and Compute Migration Stats 
#' 
#' @keywords internal
"_PACKAGE"




#
#
# 
# setwd("~/Documents/r_pack_dev/cellTracker/my_R_v4/")
# package.skeleton(name = "cellTracker", code_files = "~/Documents/r_pack_dev/cellTracker/my_R_v4/cellTracker_core.R")
# dir.create("cellTracker/data")
# source("~/Documents/r_pack_dev/cellTracker/my_R_v4/cellTracker_core.R")
# sample_tracked_cells <- load_tif(tiff_file = "~/Documents/r_pack_dev/cellTracker/input/input_sample.tif", experiment = "My Experiment 01")
# sample_tracked_cells@images$images <- sample_tracked_cells@images$images[1:10]
# sample_tracked_cells@images$attributes <- sample_tracked_cells@images$attributes[1:10]
# sample_tracked_cells@images$dim$NumberImages <- 10
# save(list = "sample_tracked_cells", file = "cellTracker/data/sample_tracked_cells.rda", compress = "xz")
# file.remove(c("cellTracker/DESCRIPTION", "cellTracker/NAMESPACE", "cellTracker/Read-and-delete-me"))
# sapply(grep("Rd$", dir("cellTracker/man/"), value = TRUE), function(x) {file.remove(paste0("cellTracker/man/", x))})
# file.copy(from = "DESCRIPTION", to = "cellTracker/")
# roxygen2::roxygenize("cellTracker/")
# file.remove(c("cellTracker/DESCRIPTION", "cellTracker/NAMESPACE"))
# file.copy(from = "NAMESPACE", to = "cellTracker/")
# file.copy(from = "DESCRIPTION", to = "cellTracker/")




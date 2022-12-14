# Class specification and constructors for the simulation class

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Validity check
valid_MizerSim <- function(object){
    errors <- character()
    validObject(object@params)
    # array dimensions
    if (length(dim(object@n)) != 3) {
	msg <- "n slot must have three dimensions"
	errors <- c(errors, msg)
    }
    if (length(dim(object@effort)) != 2) {
	msg <- "effort slot must have two dimensions"
	errors <- c(errors, msg)
    }
    if (length(dim(object@n_pp)) != 2) {
	msg <- "n_pp slot must have two dimensions"
	errors <- c(errors, msg)
    }
    # Check time dimension is good - size, dim name, and names
    if (!all(c(dim(object@n)[1],dim(object@n_pp)[1]) == dim(object@effort)[1])) {
	msg <- "First dimension of effort, n and n_pp slots must be the same length"
	errors <- c(errors, msg)
    }
    if (!all(c(names(dimnames(object@n))[1], 
               names(dimnames(object@n_pp))[1], 
               names(dimnames(object@effort))[1]) == "time")) {
	msg <- "First dimension of effort, n and n_pp slots must be called 'time'"
	errors <- c(errors, msg)
    }
    # species dimension of n
    if (dim(object@n)[2] != dim(object@params@psi)[1]) {
	msg <- "Second dimension of n slot must have same length as the species names in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n))[2] != "sp") {
	msg <- "Second dimension of n slot must be called 'sp'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@n))[2] == names(dimnames(object@params@psi))[1])) {
	msg <- "Second dimension of n slot must have same species names as in the params slot"
	errors <- c(errors, msg)
    }
    # w dimension of n
    if (dim(object@n)[3] != length(object@params@w)) {
	msg <- "Third dimension of n slot must have same length as w in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n))[3] != "w") {
	msg <- "Third dimension of n slot must be called 'w'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@n))[3] == names(dimnames(object@params@psi))[2])) {
	msg <- "Third dimension of n slot must have same size names as in the params slot"
	errors <- c(errors, msg)
    }
    # w dimension of n_pp
    if (dim(object@n_pp)[2] != length(object@params@w_full)) {
	msg <- "Second dimension of n_pp slot must have same length as w_full in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n_pp))[2] != "w") {
	msg <- "Second dimension of n_pp slot must be called 'w'"
	errors <- c(errors, msg)
    }
    if (!all(dimnames(object@n_pp)$w == names(object@params@rr_pp))) {
	msg <- "Second dimension of n_pp slot must have same size names as rr_pp in the params slot"
	errors <- c(errors, msg)
    }
    # gear dimension of effort
    if (dim(object@effort)[2] != dim(object@params@catchability)[1]) {
	msg <- "Second dimension of effort slot must have same number of gears as in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@effort))[2] != "gear") {
	msg <- "Second dimension of effort slot must be called 'gear'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@effort))[2] == names(dimnames(object@params@catchability)[1]))) {
	msg <- "Second dimension of effort slot must have same gear names as in the params slot"
	errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}

# Soundtrack: Yob - Quantum Mystic
#### Class definition ####
#' MizerSim
#' 
#' A class that holds the results of projecting a \linkS4class{MizerParams}
#' object through time.
#' 
#' \linkS4class{MizerSim} objects are created by using the \code{\link{project}} method
#' on an object of type \code{MizerParams}.
#' 
#' There are several plotting methods available to explore the contents of a
#' \code{MizerSim} object. See the package vignette for more details.
#' 
#' @slot params An object of type \linkS4class{MizerParams}.
#' @slot n Array that stores the projected community population abundances by
#'   time, species and size
#' @slot effort Array that stores the fishing effort through time by time and
#'   gear
#' @slot n_pp Array that stores the projected plankton abundance by time and
#'   size
#'   ##AAsp
#' @slot n_bb Array that stores the projected benthos abundance by time and
#'   size
#' @slot n_aa Array that stores the projected algal abundance by time and
#'   size
#' @slot diet_comp Array that stores diet composition by predator/size/prey/size 
#'   for a chosen number of steps (usually last ten years)
#'   
#' @seealso \code{\link{project}} \code{\link{MizerParams}}
#' @export

setClass(
    "MizerSim",
    representation(
        params = "MizerParams",
        n = "array",
        effort = "array",
        # AA
        n_pp = "array",
        n_bb = "array",
        n_aa = "array",
        diet_comp="array",
        temperature = "matrix",
        metTempScalar = "array",
        matTempScalar = "array",
        morTempScalar = "array",
        intTempScalar = "array",
        
        # CN adding the fleetDynamics arguments 
        effortOut = "array",
        yield = "array",
        profit = "array",
        revenue = "array",
        F = "array",
        BioOut = "list", 
        BioLimits = "list"
        
        # n_pp = "array"

    ),
    prototype = prototype(
        params = new("MizerParams"),
        n = array(
            NA,dim = c(1,1,1), dimnames = list(time = NULL, sp = NULL, w = NULL)
        ),
        effort = array(
            NA,dim = c(1,1), dimnames = list(time = NULL, gear = NULL)
        ),

        temperature = matrix(
          NA, dimnames = list(time = NULL, temperature = NULL)
        ),

        # CN again add the fleetdynamics bit
        effortOut = array(
          NA,dim = c(1,1), dimnames = list(time = NULL, gear = NULL)
        ),
        yield = array(
          NA,dim = c(1,1,1,1), dimnames = list(time = NULL, species = NULL, w = NULL, gear = NULL)
        ),
        profit = array(
          NA,dim = c(1,1), dimnames = list(time = NULL, gear = NULL)
        ),
        revenue = array(
          NA,dim = c(1,1), dimnames = list(time = NULL, gear = NULL)
        ),
        F = array(
          NA,dim = c(1,1,1,1), dimnames = list(time = NULL, species = NULL, w = NULL, gear = NULL)
        ), 
        
        BioOut = list(),
        BioLimits = list(),
        
        n_pp = array(
            NA,dim = c(1,1), dimnames = list(time = NULL, w = NULL)
        ),
        
        # AA
        n_bb = array(
          NA,dim = c(1,1), dimnames = list(time = NULL, w = NULL)
        ),
        n_aa = array(
          NA,dim = c(1,1), dimnames = list(time = NULL, w = NULL)
        ),
        diet_comp = array(
          NA, dim = c(1,1,1,1), dimnames = list( predator= NULL, pred_size = NULL, prey =NULL, prey_size=NULL)
        ),
        metTempScalar = array(
          NA, dim = c(1,1,1), dimnames = list(sp = NULL, w = NULL, temperature = NULL)
        ),
        matTempScalar = array(
          NA, dim = c(1,1,1), dimnames = list(sp = NULL, w = NULL, temperature = NULL)
        ),
        morTempScalar = array(
          NA, dim = c(1,1,1), dimnames = list(sp = NULL, w = NULL, temperature = NULL)
        ),
        intTempScalar = array(
          NA, dim = c(1,1,1), dimnames = list(sp = NULL, w = NULL, temperature = NULL)
        )
    ),
    validity = valid_MizerSim
)

setValidity("MizerSim", valid_MizerSim)
remove(valid_MizerSim)


#' Constructor for the \code{MizerSim} class
#' 
#' A constructor for the \code{MizerSim} class. This is used by the
#' \code{project} method to create \code{MizerSim} objects of the right
#' dimensions. It is not necessary for users to use this constructor.
#' 
#' @param params a \linkS4class{MizerParams} object
#' @param t_dimnames Numeric vector that is used for the time dimensions of the
#'   slots. Default = NA.
#' @param t_max The maximum time step of the simulation. Only used if t_dimnames
#'   = NA. Default value = 100.
#' @param t_save How often should the results of the simulation be stored. Only
#'   used if t_dimnames = NA. Default value = 1.
#'   
#' @return An object of type \linkS4class{MizerSim}
MizerSim <- function(params, t_dimnames = NA, t_max = 100, t_save = 1) {
  
    # If the dimnames for the time dimension not passed in, calculate them
    # from t_max and t_save
    if (any(is.na(t_dimnames))){
        t_dimnames <- seq(from = 0, to = t_max, by = t_save)
    }
    if (!is.numeric(t_dimnames)){
        stop("The t_dimnames argument must be numeric.")
    }
    if (is.unsorted(t_dimnames)) {
        stop("The t_dimnames argument should be increasing.")
    }
    no_sp <- nrow(params@species_params)
    species_names <- dimnames(params@psi)$sp
    no_w <- length(params@w)
    w_names <- dimnames(params@psi)$w
    t_dim <- length(t_dimnames)
    array_n <- array(NA, dim = c(t_dim, no_sp, no_w), 
                     dimnames = list(time = t_dimnames, 
                                     sp = species_names, w = w_names))
    
    no_gears <- dim(params@selectivity)[1]
    gear_names <- dimnames(params@selectivity)$gear
    array_effort <- array(NA, dim = c(t_dim, no_gears), 
                          dimnames = list(time = t_dimnames, 
                                          gear = gear_names))
    # temperature scalars
    matrix_temperature <- matrix(NA, nrow = t_dim, 
                          dimnames = list(time = t_dimnames, 
                                          "temperature"))
    
    # AA
    array_metTempScalar = array(NA, dim = c(dim(params@species_params)[1],length(params@w),t_dim), dimnames = list(sp = params@species_params$species, w = params@w, temperature = t_dimnames))
    array_matTempScalar = array(NA, dim = c(dim(params@species_params)[1],length(params@w),t_dim), dimnames = list(sp = params@species_params$species, w = params@w, temperature = t_dimnames))
    array_morTempScalar = array(NA, dim = c(dim(params@species_params)[1],length(params@w),t_dim), dimnames = list(sp = params@species_params$species, w = params@w, temperature = t_dimnames))
    array_intTempScalar = array(NA, dim = c(dim(params@species_params)[1],length(params@w),t_dim), dimnames = list(sp = params@species_params$species, w = params@w, temperature = t_dimnames))


    # CN add yield, profit... why so many times? 
    array_effortOut <- array(NA, dim = c(t_dim, no_gears), 
                          dimnames = list(time = t_dimnames, 
                                          gear = gear_names))
    array_yield <- array(NA, dim = c(t_dim, no_sp, no_w, no_gears), 
                         dimnames = list(time = t_dimnames,
                                         species = species_names,
                                         w = w_names,
                                         gear = gear_names))
    array_profit <- array(NA, dim = c(t_dim, no_gears), 
                          dimnames = list(time = t_dimnames, 
                                          gear = gear_names))
    array_revenue <- array(NA, dim = c(t_dim, no_gears), 
                          dimnames = list(time = t_dimnames, 
                                          gear = gear_names))
    array_F <- array(NA, dim = c(t_dim, no_sp, no_w, no_gears), 
                         dimnames = list(time = t_dimnames,
                                         species = species_names,
                                         w = w_names,
                                         gear = gear_names))
    list_BioOut <- list()
    list_BioLimits <- list()
    
    no_w_full <- length(params@w_full)
    w_full_names <- names(params@rr_pp)
    array_n_pp <- array(NA, dim = c(t_dim, no_w_full), 
                        dimnames = list(time=t_dimnames, 
                                        w = w_full_names))
    array_n_bb <- array(NA, dim = c(t_dim, no_w_full), 
                        dimnames = list(time=t_dimnames, 
                                        w = w_full_names))
    array_n_aa <- array(NA, dim = c(t_dim, no_w_full), 
                        dimnames = list(time=t_dimnames, 
                                        w = w_full_names))
    sim <- new('MizerSim',
               n = array_n, 
               effort = array_effort,
               
               # AA
               temperature = matrix_temperature,

               # CN add fleet dynamics matrices... why so many times? 
               effortOut = array_effortOut,
               yield = array_yield, 
               profit = array_profit,
               revenue = array_revenue,
               F = array_F, 
               BioOut = list_BioOut,
               BioLimits = list_BioLimits,
               
               n_pp = array_n_pp,
               n_bb = array_n_bb,
               n_aa = array_n_aa,
               params = params,
               diet_comp=as.array(1,dim = c(1,1,1,1)), #place holder for diet comp array; constructed depending on whether diet comp is requested
               metTempScalar = array_metTempScalar,
               matTempScalar = array_matTempScalar,
               morTempScalar = array_morTempScalar,
               intTempScalar = array_intTempScalar)
    return(sim)
}



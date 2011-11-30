model.specials <- function(data,specials,allowInteractions=FALSE){

  names(specials) <- specials
  Terms <- attr(data,"terms")
  spec <- lapply(specials,function(sp){
    if (length(attr(Terms,"specials")[[sp]]))
      untangle.specials(Terms,sp,1:10)
    else NULL
  })

  # -------------------------check interactions-------------------------
  
  if (allowInteractions==FALSE){
    lapply(spec[sapply(spec,length)>0],function(sp){
      ord <- attr(Terms, "order")[sp$terms]
      if (any(ord > 1))
        stop(paste(sp," can not be used in an interaction"),call.=FALSE)})
  }
  
  special.frame <- lapply(spec,function(sp){
    if (length(sp)) {
      sp.frame <- data[,sp$vars,drop=FALSE]
      names(sp.frame) <- extract.name.from.special(names(sp.frame))
      sp.frame
    }
    else NULL})
  all.varnames <- all.vars(delete.response(Terms))
  unspecified <- all.varnames[!(all.varnames %in% unlist(lapply(special.frame,names)))]
  special.frame$unspecified <- data[,unspecified,drop=FALSE]
  special.frame
}


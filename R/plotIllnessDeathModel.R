plotIllnessDeathModel <- function(stateLabels,
                                  style=1,
                                  recovery=FALSE,
                                  ...){
  if (missing(stateLabels)) labels <- c("Disease\nfree","Illness","Death")
  if (recovery==TRUE){
    idHist <- Hist(time=1:4,event=list(from=c(1,1,2,2),to=c(2,3,1,3)))
    if (style==1)
      plot(idHist,
           stateLabels=stateLabels,
           box1.row=2,
           box1.column=1,
           box2.row=1,
           box2.column=3,
           ...)
    else{
      plot(idHist,
           stateLabels=stateLabels,
           ...)
    }
  }
  else{
    idHist <- Hist(time=1:3,event=list(from=c(1,1,2),to=c(2,3,3)))
    if (style==1){
      plot(idHist,
           stateLabels=stateLabels,
           box1.row=2,
           box1.column=1,
           box2.row=1,
           box2.column=3,
           ...)
    }
    else{
      plot(idHist,
           stateLabels=stateLabels,
           ...)
    }
  }
}

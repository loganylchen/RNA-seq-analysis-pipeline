 watermarkGrob <- function(lab = "owned by Logan"){
   grob(lab=lab, cl="watermark")
 }

 ## custom draw method to
 ## calculate expansion factor on-the-fly
 drawDetails.watermark <- function(x, rot = 45, ...){
 cex <- convertUnit(unit(1,"npc"), "mm", val=TRUE) /
   convertUnit(unit(1,"grobwidth", textGrob(x$val)), "mm",val=TRUE)

 grid.text(x$lab,  rot=rot, gp=gpar(cex = cex, col="#7fcdbb",
                                        fontface = "bold", alpha = 0.5))

 }


 tidy_ensemble_id <- function(ids){
  new_ids <- c()
  for(eachid in ids){
   new_ids<-c(new_ids, strsplit(eachid,split = '.')[0])
  }
  return(new_ids)
 }
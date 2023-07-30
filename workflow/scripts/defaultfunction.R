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
  for(each_id in ids){
   new_ids<-append(new_ids, strsplit(each_id,split = '.')[1])
  }
  print(length(new_ids))
  print(length(ids))
  print(new_ids[1:min(20,length(new_ids))])
  print(ids[1:min(20,length(ids))])
  return(new_ids)
 }
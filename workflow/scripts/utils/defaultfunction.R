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
   new_ids<-append(new_ids, unlist(strsplit(each_id,split = '[.]'))[1])
  }
  print(length(new_ids))
  print(length(ids))
  print(new_ids[1:min(20,length(new_ids))])
  print(ids[1:min(20,length(ids))])
  return(new_ids)
 }


 find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
     pks <- unlist(pks)
     return(pks)
}
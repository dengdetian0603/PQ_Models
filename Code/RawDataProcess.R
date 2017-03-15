
GetTop5 <- function(rawdata, top5.names) {
  rawdata.case = subset(rawdata, Y == 1)
  rawdata.ctrl = subset(rawdata, Y == 0)
  var.names = colnames(rawdata)
  
  MBS.case = rawdata.case[, grepl("NPPCR", var.names)]
  MBS.ctrl = rawdata.ctrl[, grepl("NPPCR", var.names)]
  MSS.case = rawdata.case[, grepl("BCX", var.names)]
  
  MBS.names = substr(colnames(MBS.case), 11, 99)
  MSS.names = substr(colnames(MSS.case), 9, 99)
  
  top5.mbs.case = as.data.frame(matrix(NA, nrow = nrow(MBS.case), ncol = 5))
  top5.mbs.ctrl = as.data.frame(matrix(NA, nrow = nrow(MBS.ctrl), ncol = 5))
  top5.mss.case = as.data.frame(matrix(NA, nrow = nrow(MSS.case), ncol = 5))
  colnames(top5.mbs.case) = top5.names
  colnames(top5.mbs.ctrl) = top5.names
  colnames(top5.mss.case) = top5.names
  
  for (var in top5.names) {
    if (var %in% MSS.names) {
      top5.mss.case[, var] = MSS.case[, which(MSS.names == var)]
    } else {
      top5.mss.case[, var] = as.numeric(top5.mss.case[, var])
    }
    if (var %in% MBS.names) {
      top5.mbs.case[, var] = MBS.case[, which(MBS.names == var)]
      top5.mbs.ctrl[, var] = MBS.ctrl[, which(MBS.names == var)]
    } else {
      top5.mbs.case[, var] = as.numeric(top5.mbs.case[, var])
      top5.mbs.ctrl[, var] = as.numeric(top5.mbs.ctrl[, var])
    }
  }
  list(top5.mbs.case = top5.mbs.case,
       top5.mbs.ctrl = top5.mbs.ctrl,
       top5.mss.case = top5.mss.case)
}

GetTopK <- function(rawdata, topK.names) {
  K = length(topK.names)
  rawdata.case = subset(rawdata, Y == 1)
  rawdata.ctrl = subset(rawdata, Y == 0)
  var.names = colnames(rawdata)
  
  MBS.case = rawdata.case[, grepl("NPPCR", var.names)]
  MBS.ctrl = rawdata.ctrl[, grepl("NPPCR", var.names)]
  MSS.case = rawdata.case[, grepl("BCX", var.names)]
  
  MBS.names = substr(colnames(MBS.case), 11, 99)
  MSS.names = substr(colnames(MSS.case), 9, 99)
  
  topK.mbs.case = as.data.frame(matrix(NA, nrow = nrow(MBS.case), ncol = K))
  topK.mbs.ctrl = as.data.frame(matrix(NA, nrow = nrow(MBS.ctrl), ncol = K))
  topK.mss.case = as.data.frame(matrix(NA, nrow = nrow(MSS.case), ncol = K))
  colnames(topK.mbs.case) = topK.names
  colnames(topK.mbs.ctrl) = topK.names
  colnames(topK.mss.case) = topK.names
  
  for (var in topK.names) {
    if (var %in% MSS.names) {
      topK.mss.case[, var] = MSS.case[, which(MSS.names == var)]
    } else {
      topK.mss.case[, var] = as.numeric(topK.mss.case[, var])
    }
    if (var %in% MBS.names) {
      topK.mbs.case[, var] = MBS.case[, which(MBS.names == var)]
      topK.mbs.ctrl[, var] = MBS.ctrl[, which(MBS.names == var)]
    } else {
      topK.mbs.case[, var] = as.numeric(topK.mbs.case[, var])
      topK.mbs.ctrl[, var] = as.numeric(topK.mbs.ctrl[, var])
    }
  }
  bs.available = which(apply(topK.mbs.case, 2, function(x) mean(is.na(x))) < 0.5)
  ss.available = which(apply(topK.mss.case, 2, function(x) mean(is.na(x))) < 0.5)
  list(topK.mbs.case = topK.mbs.case,
       topK.mbs.ctrl = topK.mbs.ctrl,
       topK.mss.case = topK.mss.case,
       bs.available = bs.available,
       ss.available = ss.available)
}
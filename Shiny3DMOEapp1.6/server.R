
server<-function(input, output, session){ #code per output declared before
  
  result_auth <- secure_server(check_credentials = check_credentials(credentials))
  
  updateSelectizeInput(session=session, "gene", choices =  c(Choose = '', allAxes_detected), server = TRUE, selected="Acsm4")
  
  updateSelectizeInput(session=session, "gene1", choices =  c(Choose = '', allAxes_detected), server = TRUE, selected="Acsm4")
  
  updateSelectizeInput(session=session, "gene2", choices =  c(Choose = '', allAxes_detected), server = TRUE, selected="Nqo1")
  
  #normalize fitted log10 rpm normalized values by slice volume (assuming that the RNA is homogeneously distributed in the tissue) and IPF
  
  #multiply expression values in each slice by the fraction of tissue that corresponds to that slice
  target.LML<-reactive({
    as.numeric(fittedLML[rownames(fittedLML)==input$gene,])*sliceVolLML/(sum(sliceVolLML))
  })
  
  target.AP<-reactive({
    as.numeric(fittedAP[rownames(fittedAP)==input$gene,])*sliceVolAP/(sum(sliceVolAP))
  })
  
  target.DV<-reactive({
    as.numeric(fittedDV[rownames(fittedDV)==input$gene,])*sliceVolDV/(sum(sliceVolDV))
  })
  
  #total counts of that gene in the whole tissue (which will be an average of the values for the 3 axes)
  avGeneCounts<-reactive({
    mean(c(sum(as.numeric(logNormLML[rownames(logNormLML)==input$gene,])),  sum(as.numeric(logNormAP[rownames(logNormAP)==input$gene,])), sum(as.numeric(logNormDV[rownames(logNormDV)==input$gene,]))))
  })
  
  
  df<-reactiveValues()
  df$X<-NA
  
  observeEvent(input$topics3D,{  getTopicsMat()  })
  getTopicsMat <- function(){
    showNotification("loading...")
    df$X<-zones[,6:8]
    df$X$zone<-apply(zones[, 1:5], 1, function(x){which(x==max(x))})
  }
  
  
  M3<-reactiveValues()
  M3$X<-NA
  
  observeEvent(input$Calculate3D,{  calculate3D()  })
  calculate3D <- function(){
    showNotification("loading...")
    M3$X <- melt(ipf3D(seed=seed.3d, targetx=target.DV()*avGeneCounts()/sum(target.DV()), targety=target.LML()*avGeneCounts()/sum(target.LML()), targetz=target.AP()*avGeneCounts()/sum(target.AP()), nodecs=0))
  }
  
  #Odorants3D<-reactiveValues()
  #Odorants3D$X<-NA
  
  #observeEvent(input$searchOdorants,{  searchOdorants()  })
  #searchOdorants <- function(){
  #  showNotification("loading data...")
  #  Odorants3D$X <- read.csv("./odorants3D_dataset_ORfitVals.csv", stringsAsFactors = F, row.names=1)
  #  showNotification("ready!")
  #  output$Calculate3Dod <- renderUI({
  #    actionButton("Calculate3Dod", label = "3D pattern")
  #  })
  #}
  
  #M4<-reactiveValues()
  #M4$X<-NA
  
  #observeEvent(input$Calculate3Dod,{  calculate3Dod()  })
  #calculate3Dod <- function(){
  #  M4$X<-Odorants3D$X[which(Odorants3D$X$Odorant==input$odorant),]
  #  M4$X$value[M4$X$value<0.1]<-0
  #}
  
  ######
  
  #output$downloadData <- downloadHandler(
  #  filename = "manual.pdf",
  #  content = function(file) {
  #    file.copy("www/manual.pdf", file)
  #  }
  #)
  
  output$ThreeDStructure<-renderPlotly({
    
    palette <- colorRampPalette(c("gray", "blue", "yellow", "orange", "red"))(n=100)
    
    df<-df$X
    #df<-zones[,6:8]
    #df$zone<-apply(zones[, 1:5], 1, function(x){which(x==max(x))})
    if(!is.na(df)){
      plot_ly(df, x = ~DV, y = ~LML, z = ~AP, color = ~zone, colors = palette, marker=c(list(size=1.5))) %>% colorbar(limits = c(min(df$zone),max(df$zone)))
    }
  })
  
  output$ThreeDSlicesZones<-renderPlotly({
    
    palette <- colorRampPalette(c("gray", "blue", "yellow", "orange", "red"))(n=100)
    
    #df<-zones[,6:8]
    #df$zone<-apply(zones[, 1:5], 1, function(x){which(x==max(x))})
    df<-df$X
    
    if(!is.na(df)){
      j=1
      for(i in seq(1, max(df$AP), 3)){
        mt2<-data.frame(DV=df$DV[df$AP==i], LML=df$LML[df$AP==i], val=df$zone[df$AP==i])
        if(sum(mt2$val)==0){
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~LML, y = ~DV, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[2])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[2]*15, height=dim(seed.3d)[1]*15) %>% hide_colorbar())
        }else{
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~LML, y = ~DV, color = ~val, colors = palette, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[2])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[2]*15, height=dim(seed.3d)[1]*15) %>% colorbar(limits = c(min(df$zone),max(df$zone))) %>% hide_colorbar())
        }
        j=j+1
      }
      varlist<-ls(pattern="p[[:digit:]]")
      varlist<-mixedsort(varlist)
      plotlistLMLxDV<-list()
    
      for(i in 1:length(varlist)){
        plotlistLMLxDV[[i]]<-get(varlist[i])
      }
    
      rm(list=varlist)
      j=1
      for(i in seq(1, max(df$LML), 3)){
        mt2<-data.frame(DV=df$DV[df$LML==i], AP=df$AP[df$LML==i], val=df$zone[df$LML==i])
        if(sum(mt2$val)==0){
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~DV, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[1]*15) %>% hide_colorbar())
        }else{
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~DV, color = ~val, colors = palette, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[1]*15) %>% colorbar(limits = c(min(df$zone),max(df$zone))) %>% hide_colorbar())
        }
        j=j+1
      }
      varlist<-ls(pattern="p[[:digit:]]")
      varlist<-mixedsort(varlist)
      plotlistAPxDV<-list()
      for(i in 1:length(varlist)){
        plotlistAPxDV[[i]]<-get(varlist[i])
      }
    
      rm(list=varlist)
      j=1
      for(i in seq(1, max(df$DV), 3)){
        mt2<-data.frame(LML=df$LML[df$DV==i], AP=df$AP[df$DV==i], val=df$zone[df$DV==i])
        if(sum(mt2$val)==0){
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~LML, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(0, dim(seed.3d)[2])), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[2]*15) %>% hide_colorbar())
        }else{
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~LML, color = ~val, colors = palette, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(0, dim(seed.3d)[2])), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[2]*15) %>% colorbar(limits = c(min(df$zone),max(df$zone))) %>% hide_colorbar())
        }
        j=j+1
      }
      varlist<-ls(pattern="p[[:digit:]]")
      varlist<-mixedsort(varlist)
      plotlistAPxLML<-list()
      for(i in 1:length(varlist)){
        plotlistAPxLML[[i]]<-get(varlist[i])
      }
      projection<-switch(input$projectionZones, LMLxDV = plotlistLMLxDV, APxDV = plotlistAPxDV, APxLML = plotlistAPxLML, plotlistLMLxDV)
      if(sum(df$zone>0)){
        subplot(projection, nrows=5, titleX = T, shareX = T, shareY = T)
      }
    }
  })
  
  output$OneDPlots<-renderPlot({
    if(length(input$gene)>0){
      par(mfrow=c(2, 2))
      for(axis in c("DV", "AP", "LML")){
        if(input$gene %in% rownames(get(paste("logNorm", axis, sep="")))){
          title<-paste(input$gene, "Normalized expression across the", axis, "axis")
          target<-as.numeric(get(paste("DEGs", axis, sep=""))[rownames(get(paste("DEGs", axis, sep="")))==input$gene,])
          targetFit<-as.numeric(get(paste("fitted", axis, sep=""))[rownames(get(paste("fitted", axis, sep="")))==input$gene,])
          plot(target, xlab=paste(axis, "position"), ylab="Norm. expression", main=title)
          lines(targetFit, col="red")
        }
        else{
          plot.new()
          text(0.5, 0.2, paste(input$gene, "was not found in the", axis, "axis dataset"))
        }
      }
    }
  }, width=1000)
  
  output$download1D_DV <- downloadHandler(
    filename = function(){"DV_RPMnormExpression1D.csv"}, 
    content = function(fname){
      DVdata=logNormDV
      write.csv(DVdata, fname)
    }
  )
  
  output$download1D_AP <- downloadHandler(
    filename = function(){"AP_RPMnormExpression1D.csv"}, 
    content = function(fname){
      APdata=logNormAP
      write.csv(APdata, fname)
    }
  )
  
  output$download1D_LML <- downloadHandler(
    filename = function(){"LML_RPMnormExpression1D.csv"}, 
    content = function(fname){
      LMLdata=logNormLML
      write.csv(LMLdata, fname)
    }
  )
    
  output$ThreeDPlot<-renderPlotly({
    
    #M3=M3()
    M3=M3$X
    if(!is.na(M3)){
      names(M3)<-c("DV", "LML", "AP", "value")
      M3<-M3[-which(M$value==0),]
    
      #M3$value[M3$value<0.1]<-0
      M3$value<-log10(M3$value+1)
    
      palette1 <- colorRampPalette(c("gray", "gray", "blue", "yellow", "orange", "red"))(n = 60)
      palette2 <- colorRampPalette(c("darkred", "deeppink4", "black"))(n = 250)
      my_palette2 <- c(palette1, palette2)
      my_palette2 <- colorRampPalette(c("gray", "gray", "blue", "yellow", "orange", "red"))(n=100)
      my_palette2 <- colorRampPalette(c("gray", "blue", "yellow", "orange", "red"))(n=100)
    
    
      if(sum(M3$value>0)){
        plot_ly(M3, x = ~DV, y = ~LML, z = ~AP, color = ~value, colors = my_palette2) %>% colorbar(limits = c(min(M3$value), max(M3$value)))
      }
    }
  })
  
  output$ThreeDSlices<-renderPlotly({
    #M3=M3()
    M3=M3$X
    if(!is.na(M3)){
      names(M3)<-c("DV", "LML", "AP", "value")
      M3<-M3[-which(M$value==0),]
    
      #M3$value[M3$value<0.1]<-0
      M3$value<-log10(M3$value+1)
    
      palette1 <- colorRampPalette(c("gray", "gray", "blue", "yellow", "orange", "red"))(n = 60)
      palette2 <- colorRampPalette(c("darkred", "deeppink4", "black"))(n = 250)
      my_palette2 <- c(palette1, palette2)
      my_palette2 <- colorRampPalette(c("gray", "gray", "blue", "yellow", "orange", "red"))(n=100)
      my_palette2 <- colorRampPalette(c("gray", "blue", "yellow", "orange", "red"))(n=100)
    
      j=1
      for(i in seq(1, max(M3$AP), 3)){
        mt2<-data.frame(DV=M3$DV[M3$AP==i], LML=M3$LML[M3$AP==i], val=M3$value[M3$AP==i])
        if(sum(mt2$val)==0){
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~LML, y = ~DV, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[2])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[2]*15, height=dim(seed.3d)[1]*15) %>% hide_colorbar())
        }else{
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~LML, y = ~DV, color = ~val, colors = my_palette2, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[2])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[2]*15, height=dim(seed.3d)[1]*15) %>% colorbar(limits = c(min(M3$value), max(M3$value))) %>% hide_colorbar())
        }
        j=j+1
      }
      varlist<-ls(pattern="p[[:digit:]]")
      varlist<-mixedsort(varlist)
      plotlistLMLxDV<-list()
    
      for(i in 1:length(varlist)){
        plotlistLMLxDV[[i]]<-get(varlist[i])
      }
    
      rm(list=varlist)
      j=1
      for(i in seq(1, max(M3$LML), 3)){
        mt2<-data.frame(DV=M3$DV[M3$LML==i], AP=M3$AP[M3$LML==i], val=M3$value[M3$LML==i])
        if(sum(mt2$val)==0){
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~DV, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[1]*15) %>% hide_colorbar())
        }else{
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~DV, color = ~val, colors = my_palette2, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[1]*15) %>% colorbar(limits = c(min(M3$value), max(M3$value))) %>% hide_colorbar())
        }
        j=j+1
      }
      varlist<-ls(pattern="p[[:digit:]]")
      varlist<-mixedsort(varlist)
      plotlistAPxDV<-list()
      for(i in 1:length(varlist)){
        plotlistAPxDV[[i]]<-get(varlist[i])
      }
    
      rm(list=varlist)
      j=1
      for(i in seq(1, max(M3$DV), 3)){
        mt2<-data.frame(LML=M3$LML[M3$DV==i], AP=M3$AP[M3$DV==i], val=M3$value[M3$DV==i])
        if(sum(mt2$val)==0){
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~LML, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(0, dim(seed.3d)[2])), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[2]*15) %>% hide_colorbar())
        }else{
          assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~LML, color = ~val, colors = my_palette2, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(0, dim(seed.3d)[2])), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[2]*15) %>% colorbar(limits = c(min(M3$value), max(M3$value))) %>% hide_colorbar())
        }
        j=j+1
      }
      varlist<-ls(pattern="p[[:digit:]]")
      varlist<-mixedsort(varlist)
      plotlistAPxLML<-list()
      for(i in 1:length(varlist)){
        plotlistAPxLML[[i]]<-get(varlist[i])
      }
      projection<-switch(input$projection, LMLxDV = plotlistLMLxDV, APxDV = plotlistAPxDV, APxLML = plotlistAPxLML, plotlistLMLxDV)
      if(sum(M3$value>0)){
        subplot(projection, nrows=5, titleX = T, shareX = T, shareY = T)
      }
    }
  })
  
  output$download3D <- downloadHandler(
    filename = function(){"geneExpression3D.csv"}, 
    content = function(fname){
      M3=M3$X
      names(M3)<-c("DV", "LML", "AP", "value")
      M3<-M3[-which(M$value==0),]
      write.csv(M3, fname)
    }
  )
  
  output$geneZoneInds<-renderPlot({
  
    #par(mar=c(5, 5, 5, 20))
    barplot(as.numeric(genesZones[input$gene,]), main=input$gene, names.arg = colnames(genesZones), 
            ylab="degree of belonging", ylim=c(0, 1), cex.names = 1.2)
  
  }, width=700)
  
  output$downloadDoBs <- downloadHandler(
    filename = function(){"DoBs.csv"}, 
    content = function(fname){
      DoBs=genesZones
      write.csv(DoBs, fname)
    }
  )
  
  
  output$cellTypes<-renderPlot({
    data<-data.frame(cellType=cellTypes, log10TPMexp=log10(as.numeric(scTPMsAll[input$gene,])+1))
    data$cellType <- factor(data$cellType, levels = c("1-HBC", "2-INP1", "3-GBC", "4-mSC", "5-tHBC2", "7-iSC", "8-tHBC1", "9-iOSN", "10-INP3", "11-MVC1", "12-mOSN", "14-INP2", "15-MVC2"))
    ggplot(data, aes(x=cellType, y=log10TPMexp, fill=cellType)) + theme(axis.text.x = element_text(angle = 90)) + geom_boxplot()
  }, width=1000)
  
  output$downloadCellTypes <- downloadHandler(
    filename = function(){"scCellTypesExpressionTPM.csv"}, 
    content = function(fname){
      cellTypes=scTPMsAll
      write.csv(cellTypes, fname)
    }
  )
  
  output$twoGenesCor<-renderPlot({
    if(length(input$gene)>0){
      par(mfrow=c(2, 2))
      for(axis in c("DV", "AP", "LML")){
        if(input$gene1 %in% rownames(get(paste("logNorm", axis, sep=""))) & 
           input$gene2 %in% rownames(get(paste("logNorm", axis, sep="")))){
          title<-paste(axis, "axis")
          gene1<-as.numeric(get(paste("DEGs", axis, sep=""))[rownames(get(paste("DEGs", axis, sep="")))==input$gene1,])
          gene2<-as.numeric(get(paste("DEGs", axis, sep=""))[rownames(get(paste("DEGs", axis, sep="")))==input$gene2,])
          plotXYcorrelation(x=gene1, y=gene2, col=1, ylab=input$gene2, 
                                              xlab=input$gene1, main=title)
        }
        else{
          plot.new()
          text(0.5, 0.2, paste("A gene was not found in the", axis, "axis dataset"))
        }
      }
    }
  }, width=1000)
  
  #output$ThreeDPlotOdorants<-renderPlotly({
  #  M4=M4$X
  #  if(!is.na(M4)){
  #    my_palette2 <- colorRampPalette(c("gray", "gray", "blue", "yellow", "orange", "red"))(n = 100)
    
  #    if(sum(M4$value>0)){
  #      plot_ly(M4, x = ~DV, y = ~LML, z = ~AP, color = ~value, colors = my_palette2) %>% colorbar(limits = c(min(M4$value), max(M4$value)))
  #    }
  #  }
  #})
  
  #output$ThreeDSlicesOdorants<-renderPlotly({
    
  #  M4=M4$X
  #  if(!is.na(M4)){
    
    
  #    my_palette2 <- colorRampPalette(c("gray", "gray", "blue", "yellow", "orange", "red"))(n = 100)
    
  #    j=1
  #    for(i in seq(1, max(M4$AP), 3)){
  #      mt2<-data.frame(DV=M4$DV[M4$AP==i], LML=M4$LML[M4$AP==i], val=M4$value[M4$AP==i])
  #      if(sum(mt2$val)==0){
  #        assign(paste("p", j, sep=""), plot_ly(mt2, x = ~LML, y = ~DV, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[2])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[2]*15, height=dim(seed.3d)[1]*15) %>% hide_colorbar())
  #      }else{
  #        assign(paste("p", j, sep=""), plot_ly(mt2, x = ~LML, y = ~DV, color = ~val, colors = my_palette2, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[2])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[2]*15, height=dim(seed.3d)[1]*15) %>% colorbar(limits = c(min(M4$value), max(M4$value))) %>% hide_colorbar())
  #      }
  #      j=j+1
  #    }
  #    varlist<-ls(pattern="p[[:digit:]]")
  #    varlist<-mixedsort(varlist)
  #    plotlistLMLxDV<-list()
    
  #    for(i in 1:length(varlist)){
  #      plotlistLMLxDV[[i]]<-get(varlist[i])
  #    }
    
  #    rm(list=varlist)
  #    j=1
  #    for(i in seq(1, max(M4$LML), 3)){
  #      mt2<-data.frame(DV=M4$DV[M4$LML==i], AP=M4$AP[M4$LML==i], val=M4$value[M4$LML==i])
  #      if(sum(mt2$val)==0){
  #        assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~DV, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[1]*15) %>% hide_colorbar())
  #      }else{
  #        assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~DV, color = ~val, colors = my_palette2, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(dim(seed.3d)[1], 0)), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[1]*15) %>% colorbar(limits = c(min(M4$value), max(M4$value))) %>% hide_colorbar())
  #      }
  #      j=j+1
  #    }
  #    varlist<-ls(pattern="p[[:digit:]]")
  #    varlist<-mixedsort(varlist)
  #    plotlistAPxDV<-list()
  #    for(i in 1:length(varlist)){
  #      plotlistAPxDV[[i]]<-get(varlist[i])
  #    }
    
  #    rm(list=varlist)
  #    j=1
  #    for(i in seq(1, max(M4$DV), 3)){
  #      mt2<-data.frame(LML=M4$LML[M4$DV==i], AP=M4$AP[M4$DV==i], val=M4$value[M4$DV==i])
  #      if(sum(mt2$val)==0){
  #        assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~LML, color = 1, colors = "gray", showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(0, dim(seed.3d)[2])), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[2]*15) %>% hide_colorbar())
  #      }else{
  #        assign(paste("p", j, sep=""), plot_ly(mt2, x = ~AP, y = ~LML, color = ~val, colors = my_palette2, showlegend=FALSE) %>% layout(xaxis = list(range = c(0, dim(seed.3d)[3])), yaxis = list(range = c(0, dim(seed.3d)[2])), width=dim(seed.3d)[3]*15, height=dim(seed.3d)[2]*15) %>% colorbar(limits = c(min(M4$value), max(M4$value))) %>% hide_colorbar())
  #      }
  #      j=j+1
  #    }
  #    varlist<-ls(pattern="p[[:digit:]]")
  #    varlist<-mixedsort(varlist)
  #    plotlistAPxLML<-list()
  #    for(i in 1:length(varlist)){
  #      plotlistAPxLML[[i]]<-get(varlist[i])
  #    }
  #    projection<-switch(input$projectionOdorants, LMLxDV = plotlistLMLxDV, APxDV = plotlistAPxDV, APxLML = plotlistAPxLML, plotlistLMLxDV)
  #    if(sum(M4$value>0)){
  #      subplot(projection, nrows=5, titleX = T, shareX = T, shareY = T)
  #    }
  #  }
  #})
  
  #output$odorantZoneInds<-renderPlot({
    
  #  par(mar=c(5, 5, 5, 20))
  #  barplot(as.numeric(odorantsZones[input$odorant,]), main=input$odorant, names.arg = colnames(odorantsZones), ylab="degree of belonging", ylim=c(0, 1), cex.names = 1.2)
    
  #})
}

#### ----------------------------------------------------------------------------------- ####
##function to create clusters and calculate size, then add to linelist###
cluster_add_func <- function(linelist, contacts) {
  
  
  x <- epicontacts::make_epicontacts(linelist, contacts)
  xClust <- subset(x, cs_min = 0, cs_max = 300)
  xClust <- thin(xClust, what = "contacts")
  xClust$directed <- F
  cGraph <- cluster_fast_greedy(as.igraph(xClust))
  cGraph <- as.data.frame(cbind(cGraph$names, cGraph$membership))
  x$linelist <- merge(x$linelist,
                      cGraph,
                      by.x = "id",
                      by.y = "V1",
                      all = T)
  
  
  x$linelist$V2 <- paste0("cl_", as.character(x$linelist$V2))
  x$linelist$V2[x$linelist$V2 == "cl_NA"] <- "cl_0"
  names(x$linelist)[names(x$linelist) == 'V2'] <- 'clMembership'
  clustSize <- as.data.frame(table(x$linelist$clMembership))
  colnames(clustSize) <- c("member", "clSize")
  x$linelist <-
    merge(x$linelist, clustSize, by.x = "clMembership", by.y = "member")
  x$linelist <- x$linelist[, c(2:ncol(x$linelist), 1)]
  degs <- as.data.frame(get_degree(x))
  degs$id <- rownames(degs)
  x$linelist <- merge(x$linelist,
                      degs,
                      by.x = "id",
                      by.y = "id",
                      all = T)
  names(x$linelist)[names(x$linelist) == 'get_degree(x)'] <-
    'degrees'
  x$linelist$degrees <-
    ifelse(is.na(x$linelist$degrees), 0, x$linelist$degrees)
  x$linelist$clMembership <-
    ifelse(x$linelist$degrees == 0, "cl_NA", x$linelist$clMembership)
  x <- x[!is.na(x$linelist$clMembership)]
  lookup <- as.data.frame(unique(x$linelist$clMembership))
  lookup$cluster_number <-
    c("cl_NA", paste0("cl_", c(1:(nrow(
      lookup
    ) - 1))))
  names(lookup)[names(lookup) == 'unique(x$linelist$clMembership)'] <-
    'clMembership'
  x$linelist <-
    merge(x$linelist, lookup, by.x = "clMembership", by.y = "clMembership")
  x$linelist <- x$linelist[, -1]
  
  #rename the NA so they have a number
  x$linelist$cluster_number[x$linelist$cluster_number == "cl_NA"] = 
    paste0("cl_", length(unique(x$linelist$cluster_number)))
  
  
  return(x$linelist)
}

#### ----------------------------------------------------------------------------------- ####
### function to make tree if data is uploaded ###
plot_clusters = function(linelist, 
                         contacts, 
                         group = "id",
                         type = "network"){
  
  #make epicontacts
  x = epicontacts::make_epicontacts(linelist, contacts)
  
  if(type == "network"){
    #make epicontacts
    x1<- as.igraph(x)
    
    x1 <- simplify(x1)
    x2 <- intergraph::asNetwork(x1, amap = attrmap())
    
    p =  ggnet2(x2,
                node.color = "white",
                size="degree",
                alpha = 0.75, 
                edge.alpha = 0.5,
                legend.position  ="none")
    
    p$data<-merge(p$data,x$linelist,by.x="label",by.y="id")#ggnet is a pain, so need to add cluster degree size details manually to get a filled circle...
    p$data$size<-3+(p$data$size / max(p$data$size)) * 10##scale sizes so they are relative
    
    group_name=as.vector(p$data$cluster_number)
    id=as.vector(p$data$label)
    onset=as.vector(p$data$reported_onset_date)
    
    
    Colour = as.character(p$data[, "cluster_number"])
    
    p=p+geom_point(aes(x = p$data$x,
                       y = p$data$y,
                       fill=Colour,
                       id=id,
                       onset=onset,
                       Group=Colour), 
                   alpha = 0.7,
                   color="black",shape=21)+
      # scale_fill_manual(values = colPal)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    return(p)
    
    
  }
  
  if (type =="table" ) {
    
    knitr::kable(linelist)
    
  }
  
}


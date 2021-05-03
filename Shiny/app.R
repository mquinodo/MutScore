# 01-kmeans-app

library(MASS)
library(UniProt.ws)
library(drawProteins)
library(shiny)
library(inlmisc)
library(DT)
library(caTools)
library(shinythemes)
library(waiter)

# library(rsconnect)
# setRepositories()
# options(repos = BiocManager::repositories())
# deployApp("~/Desktop/Shiny-MutLand",appName="mutland")
# deployApp("/home/mquinodo/SYNO/MutLand/Shiny-MutLand",appName="mutscore")
# shiny::runApp("/home/mquinodo/SYNO/MutLand/Shiny-MutLand")

load(file="init.RData")

#load(file="~/Desktop/Shiny-MutLand/data-mutland.RData")


demiinterval=0.28
demigene=0.05
beforelegend=0.05
legendxfactor=1.02

ui <- navbarPage(windowTitle="MutScore",
  theme = shinytheme("cerulean"),
  title=div(img(src="IOB2.png", height="30px", width="110px"), style = "font-size:30px",strong("MutScore")),
  tabPanel(title = "Graphical representation (MutLand)",
    use_waiter(), 
    waiter_show_on_load(spin_throbber()),
    fluidRow(
      sidebarPanel(
        selectizeInput('gene', div(style = "font-size:20px", "Selected gene"), as.list(genes),selected="NSD1",options = list(maxOptions = 1000)),
        uiOutput("isoforms"),
        uiOutput("zooms"),
        radioButtons('radio', 'Type of change', choices=c('DNA', 'Protein'), selected = 'DNA', inline = TRUE),
        conditionalPanel(
          condition = "input.radio == 'DNA'",
          uiOutput("dna")
        ),
        conditionalPanel(
          condition = "input.radio == 'Protein'",
          uiOutput("protein")
        ),
        selectizeInput('constrack', 'Conservation / Score', as.list(c("MutScore","GERP++RS","PhyloP100way","PhastCons100way","pext"))),
        numericInput('conswindow', 'Moving average', 10,
                     min = 1, max = 1000),
        downloadButton('downloadPlot', 'MutLand as PDF'),
        downloadButton('downloadPlot2', HTML('MutScore distributions<br/>as PDF')),
        tags$style(type='text/css', "#downloadPlot2 { text-align: left;}"),
        width=2
      ),
      mainPanel(
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }"),
        plotOutput('plot1',height="850px",width="1180px"),
        plotOutput('plot2',height="600px",width="1180px")
      )
    )
  ),
  tabPanel(title = "Missense variants",
     fluidRow(
       titlePanel(h2(textOutput('titlevariant'),align="center"))
     ),
     br(),
    fluidRow(
      sidebarPanel(
        checkboxGroupInput("types","Type of variants:",
                           c("Pathogenic and likely pathogenic (PLP)" = "PLP",
                             "Benign and liekly benign (BLB)" = "BLB",
                             "Variants of unknown significance (VUS)" = "VUS",
                             "Conflicting interpretation (CI)" = "CON",
                             "gnomAD frequent variants (gnomAD)" = "gnomAD"
                             ),
                           selected=c("PLP")
        ),
        downloadButton('download1', 'Download data as tsv'),
        width=2
      ),
      column(9,
        DT::dataTableOutput('table1')
      )
    )
  ),
  tabPanel(title = "Clusters",
   fluidRow(
     titlePanel(h2(textOutput('titlecluster'),align="center"))
   ),
   titlePanel(h3("Regions with enriched PLP variants",align="center")),
   fluidRow(
     column(10,offset=1,
      DT::dataTableOutput('table2')
     ),
   ),
   titlePanel(h3("Regions with enriched BLB variants",align="center")),
   fluidRow(
     column(10,offset=1,
      DT::dataTableOutput('table3')
     ),
   ),
   titlePanel(h3("Overall PLP clustering",align="center")),
   fluidRow(
     column(10,offset=1,
      DT::dataTableOutput('table4')
     ),
   ),
   titlePanel(h3("Overall BLB clustering",align="center")),
   fluidRow(
     column(10,offset=1,
      DT::dataTableOutput('table5')
     ),
   )
   ),
  tabPanel(title = "Overall clustering",
   fluidRow(
    sidebarPanel(
      div(style = "font-size:20px", "Downloads"),
      br(),
      p('PLP clustering data'),
        downloadButton('download2', 'Download data as tsv'),
      br(),
      br(),
      p('BLB clustering data'),
      downloadButton('download3', 'Download data as tsv'),
        width=2
    ),
    mainPanel(
      h3("Clustering of PLP variants per isoform",align="center"),
      DT::dataTableOutput('table6'),
      br(),
      h3("Clustering of BLB variants per isoform",align="center"),
      DT::dataTableOutput('table7')
    )
   )
  ),
  
  tabPanel(title = "About - contact",
           titlePanel(h3("About")),
           p("Large-scale sequencing efforts of genomes from patients often lead to the identification of DNA variants of unknown significance (VUS), i.e. of genomic variations that cannot be immediately recognized as pathogenic mutations or benign DNA changes. Despite many in silico predictors exist, none of them takes advantage from the wealth of information contained in the positional clustering of mutations already detected in disease-associated genes."),
             p("We have therefore developed MutScore, a new pathogenicity predictor that integrates unsupervised features of single DNA variants with information derived from such clusters. The predictive model for MutScore was trained with a random forest approach on medically-relevant mutations and subsequently tested against various genomic databases for both hereditary conditions and cancer (ClinVar, HGMD, and DoCM), achieving a very high performance."),
             p("The use of MutScore on 840 genes from the ClinVar database also allowed the detection of significant clusters of disease-associated and of benign variants in 505 and 345 of them, respectively, revealing protein domains with diverging functional importance. In addition, an open-access web-based application, MutLand, was developed to provide a comprehensive graphical landscape of all known medically-relevant and clearly benign DNA variants for individual genes, as a help in appraising new VUS identified in clinical testing."),
             p("Altogether, our work reveals the presence of widespread clustering of missense variants associated with normal and clinical phenotypes and that this information can be systematically used to improve and to understand pathogenicity at the molecular level."),
           br(),
           titlePanel(h3("Contact")),
           p("Mathieu Quinodoz"),
           div("mathieu.quinodoz[at]iob.ch", style = "color:blue"),
           br(),
           p("Institute of Molecular and Clinical Ophthalmology Basel (IOB)"),
           p("Mittlere Strasse 91, 4031 Basel, Switzerland"),
           a("IOB website", href="https://www.iob.ch/"),
           br(),
           img(src = "IOB2.png", height="90px", width="330px")
  ),
  
  tabPanel(title = "Download",
           titlePanel(h3("MutScore")),
           br(),
           p("MutScore and MutLand are freely available for non-commercial use. For other uses, please contact IOB."),
           br(),
           p("MutScore values for all missense variants for hg19 genome build can be downloaded on our github page:"),
           a("github.com/mquinodo/MutScore", href="https://github.com/mquinodo/MutScore", target="_blank")
  ),
  
  tabPanel(title = "How to cite",
           titlePanel(h3("Manuscript")),
           br(),
           p("MutScore / MutLand manuscript will be available soon."),
  ),
  
  tabPanel(title = div(strong("MutScore-batch"),style = "font-size:20px"),
           titlePanel(h3("MutScore-batch")),
           br(),
           p("For MutScore-batch, please visit this webpage:"),
           a("iob-genetic.shinyapps.io/mutscore-batch", href="https://iob-genetic.shinyapps.io/mutscore-batch", target="_blank")
  )
)

server <- function(input, output, session) {

  print("Loading initial databases")
  time1=Sys.time()
  load(file="data-mutland-small.RData")
  load(file="cut/isosel.RData")
  load(file="maxiso-all.RData")
  load(file=paste("cut/data-111.RData",sep=""))
  data=dataallgenes
  scores=scores[,c(1,2,3,4,5,6,8,7),]
  waiter_hide()
  seldna=""
  selprot=""
  globaltemp=""
  time2=Sys.time()
  print(time2-time1)
  
  
  variables = reactiveValues(mut2=mut1,dataallgenes2=dataallgenes,scores2=scores,data12=data1,isoDNA="",isoPROT="")
  
  vals <- reactiveValues()
  
  output$isoforms <- renderUI({
    print("Changing output isoform")
    if(input$gene!=""){
      print("Changing output isoform done")
      isos=geneinfo[which(geneinfo[,2]==input$gene & geneinfo[,1]!="."),1]
      print(as.list(names(sort(table(isos),decreasing=TRUE))))
      selectInput("isoform", "Selected isoform", as.list(names(sort(table(isos),decreasing=TRUE))))
    }
  })
  
   observeEvent(input$isoform, {
     print("Observe isoform and loading mutscore")
     time1=Sys.time()
     if(length(input$isoform)>0 & input$isoform!=""){
        i=as.numeric(isosel[which(isosel[,1]==input$isoform),2])
        load(file=paste("cut/data-",i,".RData",sep=""))
        data=dataallgenes

        scores=scores[,c(1,2,3,4,5,6,8,7),]
        variables$mut2=unique(mut1[which(mut1[,1]==input$isoform),])
        variables$scores2=scores
        variables$data12=data1
        variables$dataallgenes2=dataallgenes

        time2=Sys.time()
        print(time2-time1)
     }
     else {print("Not loaded")}
  })

  output$zooms <- renderUI({
    maxpos=as.numeric(maxiso[which(maxiso[,1]==input$isoform),2])
    if(length(maxpos)==0){maxpos=1}
    sliderInput("zoom", "Zoom range", min=0, max=maxpos, value=c(0,maxpos))
  })
  
  output$dna <- renderUI({
    print("Updating DNA")
    
    bla=dataInput2()
    isoform=bla[1]
    
    print(paste("Status :",bla[1]," ",bla[2]," ",bla[3]," ",bla[4]," ",sep=""))
    
    mut1=variables$mut2
    cpos1=mut1[which(mut1[,1]==isoform),2]
    cpos1=cpos1[order(mut1[which(mut1[,1]==isoform),3])]
    cpos1=c("-",cpos1)
    
    seldna=bla[2]
    
    if(length(bla)>3){
      print("Yes")
      if(bla[4]=="Protein"){
        print("YES")
        seldna="-"
      }
    }

    selectInput("cpos", "DNA change (red line)", as.list(cpos1),selected=seldna)
  })
  
  output$protein <- renderUI({
    print("Updating protein")
    mut1=variables$mut2
    
    bla=dataInput2()
    isoform=bla[1]
      
    ppos1=mut1[which(mut1[,1]==isoform ),4]
    ppos1=ppos1[order(mut1[which(mut1[,1]==isoform),3])]
    ppos1=c("-",ppos1)
    ppos1=unique(ppos1)
    
    selprot=input$ppos
    
    if(length(bla)>3){
      if(bla[4]=="DNA"){
        selprot="-"
      }
    }
    
    if(variables$isoPROT!=isoform){selprot="-"}
    variables$isoPROT=isoform
    
    selectInput("ppos", "Protein change (red line)", as.list(ppos1),selected=selprot)
  })
  
  
  output$table1 <- DT::renderDataTable({
    data1sel=data1[which(data1[,2]==input$isoform),]
    data1sel=data1sel[which(is.element(data1sel[,3],input$types)),]
    a=sort(data1sel[,4],index.return=TRUE)$ix
    data1sel=data1sel[a,]
    data1sel[which(data1sel[,3]=="CON"),3]="CI"
    colnames(data1sel)=c("Gene","Isoform","Category","Amino acid position","Protein change")
    DT::datatable(data1sel,options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 25,columnDefs = list(list(className = 'dt-center', targets =2:3))), rownames=FALSE)
  })

  output$download1 <- downloadHandler(
    filename = function() { paste("Missense-variants_",input$gene,"_",input$isoform,'.tsv', sep='') },
    content = function(file) {
      data1sel=data1[which(data1[,2]==input$isoform),]
      data1sel=data1sel[which(is.element(data1sel[,3],input$types)),]
      a=sort(data1sel[,4],index.return=TRUE)$ix
      data1sel=data1sel[a,]
      data1sel[which(data1sel[,3]=="CON"),3]="CI"
      colnames(data1sel)=c("Gene","Isoform","Category","Amino acid position","Protein change")
      write.table(data1sel, file, row.names = FALSE,sep="\t",quote=FALSE)
    }
  )
  
  output$table2 <- DT::renderDataTable({
    colnames(regionsp1)=c("Gene","Isoform","Begin","End","Nb PLP inside","Nb BLB inside")
    reg=regionsp1[which(regionsp1[,2]==input$isoform),]
    DT::datatable(reg, options = list(paging = FALSE,searching=FALSE,columnDefs = list(list(className = 'dt-center', targets =2:5))))
  })
  
  output$table3 <- DT::renderDataTable({
    colnames(regionsb1)=c("Gene","Isoform","Begin","End","Nb PLP inside","Nb BLB inside")
    reg3=regionsb1[which(regionsb1[,2]==input$isoform),]
    DT::datatable(reg3, options = list(paging = FALSE,searching=FALSE,columnDefs = list(list(className = 'dt-center', targets =2:5))))
  })
  
  output$table4 <- DT::renderDataTable({
    colnames(resip1)=c("Gene","Isoform","Part of PLPs inside","Part of regions in gene","Clustering score","Odds-ratio","P-value (Fisher)")
    reg4=resip1[which(resip1[,2]==input$isoform),]
    if(is.null(dim(reg4))){reg4=t(as.matrix(reg4))}
    DT::datatable(reg4, options = list(paging = FALSE,searching=FALSE,columnDefs = list(list(className = 'dt-center', targets =2:6))))
  })
  
  output$table5 <- DT::renderDataTable({
    colnames(resib1)=c("Gene","Isoform","Part of BLBs inside","Part of regions in gene","Clustering score","Odds-ratio","P-value (Fisher)")
    reg5=resib1[which(resib1[,2]==input$isoform),]
    if(is.null(dim(reg5))){reg5=t(as.matrix(reg5))}
    DT::datatable(reg5, options = list(paging = FALSE,searching=FALSE,columnDefs = list(list(className = 'dt-center', targets =2:6))))
  })
  
  output$table6 <- DT::renderDataTable({
    colnames(resi1)=c("Gene","Isoform","Part of PLPs inside","Part of regions in gene","Clustering score","Odds-ratio","P-value (Fisher)","FDR corrected p-value")
    resi1=resi1[order(resi1[,5],decreasing=TRUE),]
    resi1=as.data.frame(resi1)
    resi1[,7]=as.numeric(resi1[,7])
    resi1[,8]=as.numeric(resi1[,8])
    DT::datatable(resi1,options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10,columnDefs = list(list(className = 'dt-center', targets =2:7))), rownames=FALSE)
  })
  
  output$table7 <- DT::renderDataTable({
    colnames(resi2)=c("Gene","Isoform","Part of BLBs inside","Part of regions in gene","Clustering score","Odds-ratio","P-value (Fisher)","FDR corrected p-value")
    resi2=resi2[order(resi2[,5],decreasing=TRUE),]
    resi2=as.data.frame(resi2)
    resi2[,7]=as.numeric(resi2[,7])
    resi2[,8]=as.numeric(resi2[,8])
    DT::datatable(resi2,options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10,columnDefs = list(list(className = 'dt-center', targets =2:7))), rownames=FALSE)
  })
  
  output$download2 <- downloadHandler(
    filename = function() { 'PLP-clustering.tsv' },
    content = function(file) {
      colnames(resi1)=c("Gene","Isoform","Part of PLPs inside","Part of regions in gene","Clustering score","Odds-ratio","P-value (Fisher)","FDR corrected p-value")
      resi1=resi1[order(resi1[,5],decreasing=TRUE),]
      resi1=as.data.frame(resi1)
      resi1[,7]=as.numeric(resi1[,7])
      resi1[,8]=as.numeric(resi1[,8])
     write.table(resi1, file, row.names = FALSE,sep="\t",quote=FALSE)
    }
  )
  
  output$download3 <- downloadHandler(
    filename = function() { 'BLB-clustering.tsv' },
    content = function(file) {
      colnames(resi2)=c("Gene","Isoform","Part of BLBs inside","Part of regions in gene","Clustering score","Odds-ratio","P-value (Fisher)","FDR corrected p-value")
      resi2=resi2[order(resi2[,5],decreasing=TRUE),]
      resi2=as.data.frame(resi2)
      resi2[,7]=as.numeric(resi2[,7])
      resi2[,8]=as.numeric(resi2[,8])
      write.table(resi2, file, row.names = FALSE,sep="\t",quote=FALSE)
    }
  )
  
  output$titlecluster <- renderText({
    paste("Clustering information for ",input$gene," (",input$isoform,")",sep="")
  })
  
  output$titlevariant <- renderText({
    paste("Missense variants for ",input$gene," (",input$isoform,")",sep="")
  })
  
  dataInput <- reactive({
    print("Reactive datainput")
    time1=Sys.time()
    c(isolate(input$isoform),input$zoom,input$constrack,input$conswindow,input$cpos,input$ppos,isolate(input$radio))
  })
  
  dataInput2 <- reactive({
    print("Reactive datainput2")
    time1=Sys.time()
    c(input$isoform,input$cpos,input$ppos,input$radio)
  })

  plotInput <- reactive({
    
    print("DOING PLOT!!!!")
    time1=Sys.time()
    mut1=variables$mut2
    data1=variables$data12
    dataallgenes=variables$dataallgenes2
    scores=variables$scores2
    
    bla=dataInput()
    
    isoform=bla[1]
    zoom=c(as.numeric(bla[2]),as.numeric(bla[3]))
    constrack=bla[4]
    conswindow=as.numeric(bla[5])
    cchange=bla[6]
    pchange=bla[7]
    radio=bla[8]
    
    data=dataallgenes[which(dataallgenes[,2]==isoform),]
    
    print(paste("Check plot: ",pchange,", ",cchange,", ",radio,", ",dim(data)[1],sep=""))
    print(isoform)

    if(is.na(pchange)==FALSE & is.na(cchange)==FALSE & is.na(radio)==FALSE & dim(data)[1]>0 & zoom[2]>1){
    
    if(is.na(pchange)){pchange="-"}
    if(is.na(cchange)){cchange="-"}
    if(is.na(radio)){radio="DNA"}
    
    if(pchange!="-" & radio=="Protein"){
      cchange=mut1[which(mut1[,1]==isoform & mut1[,4]==pchange)[1],2]
    }
    
    if(pchange=="-" & radio=="Protein"){
      cchange="-"
    }
    
    if(cchange!="-" & radio=="DNA"){
      pchange=mut1[which(mut1[,1]==isoform & mut1[,2]==cchange),4]
    }
    
    if(cchange=="-" & radio=="DNA"){
      pchange="-"
    }

    maximum=as.numeric(maxiso[which(maxiso[,1]==isoform),2])
    
    data=dataallgenes[which(dataallgenes[,2]==isoform),]
    
    gene=as.character(data[1,5])
    
    # other types of variants
    datanc1=dataallgenes[which(dataallgenes[,5]==gene),]
    datanc=datanc1[which(datanc1[,4]=="intronic" | datanc1[,4]=="UTR3" | datanc1[,4]=="UTR5" | datanc1[,4]=="splicing"),]
    data=data[which(data[,4]=="frameshift deletion" | data[,4]=="frameshift insertion" | data[,4]=="stopgain" | data[,4]=="nonsynonymous SNV" | data[,4]=="synonymous SNV"),]
    
    data=data[which(data$pos>=zoom[1] & data$pos<=zoom[2]),]
    data$pos=round((data$pos-zoom[1])/(zoom[2]-zoom[1])*maximum)
    
    # choosing frequent gnomAD variants for ClinVar plot
    data[which(data[,14]==0),14]=(-6)
    data[,14]=as.numeric(data[,14])
    data[which(data[,15]=="-7"),1]="gnohigh"
    
    # separating missenses and LoF
    dataall=data
    data2=data[which(data[,4]=="frameshift deletion" | data[,4]=="frameshift insertion" | data[,4]=="stopgain"),]
    data=data[which(data[,4]=="nonsynonymous SNV"),]
    databoth=rbind(data,data2)

    # computing steps for graphics of clinvar and gnomad
    {
      # for plot 1
      
      PLP=which(data[,1]=="PLP")
      BLB=which(data[,1]=="BLB")
      gnohigh=which(data[,1]=="gnohigh")
      m=max(table(data[PLP,18]),table(data[BLB,18]),table(data[gnohigh,18]))
      step1=demiinterval/m
      m1=m
      if(m1>100000){m1=1}
      if(m1<1){m1=1}
      
      # for plot 2
      PLP=which(data2[,1]=="PLP")
      BLB=which(data2[,1]=="BLB")
      gnohigh=which(data2[,1]=="gnohigh")
      m=0
      if (dim(data2)[1]>0){
        m=max(table(data2[PLP,18]),table(data2[BLB,18]),table(data2[gnohigh,18]))
      } else {
        m=1
      }
      step2=demiinterval/m
      m2=m
      step2=min(step1,step2)
      m2=max(m1,m2)
      if(m2>100000){m2=1}
      if(m2<1){m2=1}
      
      # for plot 3
      #conservation
      nb=3
      if(constrack=="GERP++RS"){nb=3}
      if(constrack=="PhyloP100way"){nb=4}
      if(constrack=="PhastCons100way"){nb=5}
      if(constrack=="pext"){nb=6}
      if(constrack=="MutScore"){nb=8}

      size=max(maximum,1)
      cons=scores[which(scores[,1]==isoform),]
      
      a0=as.matrix(cons[,c(7,..nb)])
      a1=1:size
      a1[1:size]=NA
      for (i in 1:size){
        a1[i]=mean(a0[which(a0[,1]==i),2])
      }
      
      consw=runmean(a1,conswindow)
        #window(as.matrix(cons[,c(7,..nb)]),start=1,end=size,deltat=conswindow)

      ma=ceiling(max(consw,na.rm=T))
      mi=floor(min(0,min(consw,na.rm=T)))
      if(mi==ma){mi=mi-1}
      m3=ma-mi
      step3=demiinterval/m3
      fact=1/step3
      consw=3.5+demigene+(consw-mi)/fact
      
      # for plot 4
      VUS=which(data[,1]=="VUS")
      CON=which(data[,1]=="CON")
      m=0
      if (length(VUS)+length(CON)>0){
        m=max(table(data[VUS,18]),table(data[CON,18]))
      } else {
        m=1
      }
      mtemp=m
      m=0
      VUS=which(data2[,1]=="VUS")
      CON=which(data2[,1]=="CON")
      
      if (dim(data2)[1]>0){
        m=max(table(data2[VUS,18]),table(data2[CON,18]))
      } else {
        m=1
      }
      m4=max(m,mtemp)
      step4=demiinterval/m4
      
      if(m3>100000){m3=1}
      if(m3<1){m3=1}
      
      if(m4>100000){m4=1}
      if(m4<1){m4=1}
      
    }
    print(m1)
    print(m2)
    print(m3)
    print(m4)
    
    defaultW <- getOption("warn") 
    options(warn = -1) 
    nm=strsplit(as.character(isoform),"\\.")[[1]][1]
    options(warn = defaultW)
    enst=conv[which(conv[,4]==nm),2][1]
    
    #layout(matrix(c(1,1,1,1,2,3), nrow = 6, ncol = 1, byrow = TRUE))
    # begin to plot
    plot(0,1000,ylim=c(2,6),col=2,xlim=c(0,max(as.numeric(maxiso[which(maxiso[,1]==isoform),2]),1)*1.32), xlab="",ylab="",
         main=bquote(bold("MutLand plot for " ~ bolditalic(.(gene)) ~ " (" ~ .(as.character(isoform)) ~ " - " ~ .(as.character(enst)) ~ ")")),yaxt='n',xaxt='n',cex.main=2,yaxs="i",cex.lab=2.5,font.main=2)
    mtext("Amino acid position",side=1,line=3,cex=1.5,at=max(maximum,1)/2,adj=0.5)
    
    if(cchange!="-"){
      position=mut1[which(mut1[,1]==isoform & mut1[,2]==cchange),5]
      posiline=round((position-zoom[1])/(zoom[2]-zoom[1])*maximum)
      if (position>0){abline(v=posiline,col='red')}
    }
    
    #Grids and rectangle
    {
      # vertical grid
      # choosing size between x-axis labels
      po=c(10,50,100,200,500,1000,2000,4000,6000)
      a=po[which.min(abs(max(zoom[2]-zoom[1],1)/po-10))]
      b=0:20
      ax2=a*b+round(zoom[1]/a)*a
      
      ax=(ax2-zoom[1])/(zoom[2]-zoom[1])*maximum
      
      rem=which(ax2<=zoom[2] & ax2>=zoom[1])
      ax2=ax2[rem]
      ax=ax[rem]
      
      print("Debug")
      print(zoom[1])
      print(zoom[2])
      print(maximum)
      
      axis(1,at=ax,labels=ax2,cex.axis=1.5)
      for (i in 1:length(ax)){
        abline(v=ax[i],lty=3,col='lightgray')
      }
      
      # horizontal grids
      y=c()
      lab=c()
      for (i in 1:m1){
        y=c(y,5.5-demigene-step1*i,5.5+demigene+step1*i)
        lab=c(lab,i,i)
      }
      for (i in 1:m2){
        y=c(y,4.5-demigene-step2*i,4.5+demigene+step2*i)
        lab=c(lab,i,i)
      }
      for (i in 0:m3){
        y=c(y,3.5+demigene+step3*i)
        lab=c(lab,i+mi)
      }
      for (i in 1:m4){
        y=c(y,2.5-demigene-step4*i,2.5+demigene+step4*i)
        lab=c(lab,i,i)
      }
      axis(2, at = y,lab,cex.axis=0.6)
      for (i in 1:length(y)){
        abline(h=y[i],lty=3,col='lightgray')
      }
      abline(h=1)
      abline(h=2)
      abline(h=3)
      abline(h=4)
      abline(h=5)
      
      # gene rectangles
      rect(0,5.5-demigene,max(maximum,1),5.5+demigene,col='lightgray')
      rect(0,4.5-demigene,max(maximum,1),4.5+demigene,col='lightgray')
      rect(0,3.5-0.03,max(maximum,1),3.5+0.03,col='lightgray')
      rect(0,2.5-demigene,max(maximum,1),2.5+demigene,col='lightgray')
      rect(0,1.5-demigene,max(maximum,1),1.5+demigene,col='lightgray')
      rect(0,0.5-demigene,max(maximum,1),0.5+demigene,col='lightgray')
      
    }
    
    # legends
    {
      # plot 1
      text(max(maximum,1)*legendxfactor*(-0.03),6-beforelegend*1.5,"ClinVar missense variants",adj=0,cex=1.5,font=2)
      legend(max(maximum,1)*legendxfactor,6-beforelegend*1,
             c("ClinVar* missense PLP variants","ClinVar* missense BLB variants","Frequent** gnomAD missense variants "),
             col=c('orange',3,4),pch=c(19,19,19),lty=NA,lwd=4,cex=1,bg='white', pt.cex=1.2)
      text(max(maximum,1)*legendxfactor*1.02,2.1,"** with AF > max(AF of PLPs)]",adj=0,cex=1.2)
      text(max(maximum,1)*legendxfactor*1.02,2.2, "* ClinVar version of 21.11.2020",adj=0,cex=1.2)
      mtext("Missense variants per aa",side=2,line=3,srt=90,cex=1.2,at=5.5)
      nb1=length(which(data[,1]=="PLP"))
      text(max(maximum,1)/2,6-beforelegend*1.5,paste("PLP missense variants from ClinVar (N=",nb1,")",sep=""),adj=0.5,cex=1.2)
      nb1=length(which(data[,1]=="BLB"))
      text(max(maximum,1)/2,5+beforelegend*1.5,paste("BLB missense variants from ClinVar (N=",nb1,")",sep=""),adj=0.5,cex=1.2)
      
      nb1=length(which(datanc[,4]=="intronic" & datanc[,1]=="PLP"))
      nb2=length(which(datanc[,4]=="UTR5" & datanc[,1]=="PLP"))
      nb3=length(which(datanc[,4]=="UTR3" & datanc[,1]=="PLP"))
      nb4=length(which(dataall[,4]=="synonymous SNV" & dataall[,1]=="PLP"))
      if(nb1+nb2+nb3+nb4>0){
        text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*9,"In ClinVar there are also:",adj=0,cex=1.2)
        bot=11.5
        if(nb1>0){
          text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*bot,paste("- ",nb1, " intronic PLP variant(s)",sep=""),adj=0,cex=1)
          bot=bot+2
        }
        if(nb2>0){
          text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*bot,paste("- ",nb2, " 5'UTR PLP variant(s)",sep=""),adj=0,cex=1)
          bot=bot+2
        }
        if(nb3>0){
          text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*bot,paste("- ",nb3, " 3'UTR PLP variant(s)",sep=""),adj=0,cex=1)
          bot=bot+2
        }
        if(nb4>0){
          text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*bot,paste("- ",nb4, " synonymous PLP variant(s)",sep=""),adj=0,cex=1)
          bot=bot+2
        }
      } else {
        text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*9,"No non-coding or synonymous PLP",adj=0,cex=1)
      }
      
      # exons
      
      exons=unique(cons[,c(2,7)])
      nbexons=max(exons[,1])
      minexon=min(exons[,1])
      posi3=1
      text(0,5.5,"Exons ",adj=1,cex=0.8)
      for (i in minexon:nbexons){
        
        if(length(which(exons[,1]==i))>0){

          posi=max(exons[which(exons[,1]==i),2])
          
          posi2=round((posi-zoom[1])/(zoom[2]-zoom[1])*maximum)
          if(posi2>=0 & posi2<=maximum){
            rect(posi2,5.5-demigene,posi2,5.5+demigene)
          }
          if((posi2-posi3)>max(maximum,1)/50 & (posi2-posi3)/2+posi3>maximum*0.01 & (posi2-posi3)/2+posi3<maximum*0.99){
            text((posi2-posi3)/2+posi3,5.5,i,adj=0.5,cex=0.8)
          }
          if(posi2=="-Inf"){posi2=1}
          posi3=posi2
        
        }
      }
      
 
      # plot 2
      text(max(maximum,1)*legendxfactor*(-0.03),5-beforelegend*1.5,"Clinvar LoF variants",adj=0,cex=1.5,font=2)
      legend(max(maximum,1)*legendxfactor,5-beforelegend*1,
             c("ClinVar* LoF variants","ClinVar* LoF BLB variants","Frequent** gnomAD LoF variants"),
             col=c(2,'darkgreen','darkviolet'),pch=c(19,19,19),lty=NA,lwd=4,cex=1,bg='white',pt.cex=1.2)
      mtext("LoF variants per aa",side=2,line=3,srt=90,cex=1.2,at=4.5)
      nb2=length(which(data2[,1]=="PLP"))
      text(max(maximum,1)/2,5-beforelegend*1.5,paste("PLP LoF variants from ClinVar (N=",nb2,")",sep=""),adj=0.5,cex=1.2)
      nb2=length(which(data2[,1]=="BLB"))
      text(max(maximum,1)/2,4+beforelegend*1.5,paste("BLB LoF variants from ClinVar (N=",nb2,")",sep=""),adj=0.5,cex=1.2)
      
      nb1=length(which(datanc[,4]=="splicing" & datanc[,1]=="PLP"))
      if (nb1>0){
        text(max(maximum,1)*legendxfactor*1.02,5-beforelegend*9,"In ClinVar there are also:",adj=0,cex=1.2)
        text(max(maximum,1)*legendxfactor*1.02,5-beforelegend*11.5,paste("- ",nb1, " splicing PLP variant(s)",sep=""),adj=0,cex=1)
      } else {
        text(max(maximum,1)*legendxfactor*1.02,5-beforelegend*9,"No splicing PLP reported",adj=0,cex=1)
      }
      
      # exons
      
      exons=unique(cons[,c(2,7)])
      nbexons=max(exons[,1])
      minexon=min(exons[,1])
      posi3=1
      text(0,4.5,"Exons ",adj=1,cex=0.8)
      for (i in minexon:nbexons){
        if(length(which(exons[,1]==i))>0){
          posi=max(exons[which(exons[,1]==i),2])
          posi2=round((posi-zoom[1])/(zoom[2]-zoom[1])*maximum)
          if(posi2>=0 & posi2<=maximum){
            rect(posi2,4.5-demigene,posi2,4.5+demigene)
          }
          if((posi2-posi3)>max(maximum,1)/50 & (posi2-posi3)/2+posi3>maximum*0.01 & (posi2-posi3)/2+posi3<maximum*0.99){
            text((posi2-posi3)/2+posi3,4.5,i,adj=0.5,cex=0.8)
          }
          if(posi2=="-Inf"){posi2=1}
          posi3=posi2
        }
      }
      
      # plot 3
      legend(max(maximum,1)*legendxfactor,4-beforelegend*2,paste(constrack," ","***",sep=""),col=3,lty=1,lwd=4,cex=1,bg='white')
      mtext("Score",side=2,line=3,srt=90,cex=1.2,at=3.65)
      text(max(maximum,1)*legendxfactor*(-0.03),4-beforelegend*1.5,"Protein domains and scores",cex=1.5,font=2,adj=0)
      text(max(maximum,1)*legendxfactor*1.02,4-beforelegend*6.5, paste("*** averaged with a window of size ",conswindow,sep=""),adj=0,cex=0.8)
      
      if(length(which(data[,1]=="PLP"))<10){
        text(max(maximum,1)*legendxfactor*1.02,4-beforelegend*8,"Note: MutScore computed without positional score",adj=0,cex=0.8)
      }
      
      if(cchange!="-"){
        text(max(maximum,1)*legendxfactor,4-beforelegend*13, paste("",cchange," - ",pchange,sep=""),adj=0,cex=1.7,col=2,font=2)
        muto=mut1[which(mut1[,2]==cchange & mut1[,1]==isoform),6]
        text(max(maximum,1)*legendxfactor,4-beforelegend*16, paste("MutScore = ",muto,sep=""),adj=0,cex=1.7,col=2,font=2)
        rect(max(maximum,1)*legendxfactor*0.99, 4-beforelegend*18, max(maximum,1)*legendxfactor*1.3, 4-beforelegend*11,border=2,lwd=2)
      }
      
      # exons
      
      exons=unique(cons[,c(2,7)])
      nbexons=max(exons[,1])
      minexon=min(exons[,1])
      posi3=1
      text(0,2.5,"Exons ",adj=1,cex=0.8)
      for (i in minexon:nbexons){
        if(length(which(exons[,1]==i))>0){
          posi=max(exons[which(exons[,1]==i),2])
          posi2=round((posi-zoom[1])/(zoom[2]-zoom[1])*maximum)
          if(posi2>=0 & posi2<=maximum){
            rect(posi2,2.5-demigene,posi2,2.5+demigene)
          }
          if((posi2-posi3)>max(maximum,1)/50 & (posi2-posi3)/2+posi3>maximum*0.01 & (posi2-posi3)/2+posi3<maximum*0.99){
            text((posi2-posi3)/2+posi3,2.5,i,adj=0.5,cex=0.8)
          }
          if(posi2=="-Inf"){posi2=1}
          posi3=posi2
        }
      }
      
      # plot 6
      text(max(maximum,1)*legendxfactor*(-0.03),3-beforelegend*1.5,"VUS and CI variants",adj=0,cex=1.5,font=2)
      label=c("ClinVar* VUS missense","ClinVar* CI missense","ClinVar* VUS LoF","ClinVar* CI LoF")
      text(max(maximum,1)*legendxfactor*1.02,3-beforelegend*8.2, "* ClinVar version of 21.11.2020",adj=0,cex=0.8)
      legend(max(maximum,1)*legendxfactor,3-beforelegend*4,label,col=c(4,3,'orange',2),lty=NA,lwd=4,pch=c(19,19,19,19),cex=1,bg='white',pt.cex=1.2)
      mtext("Variants per aa",side=2,line=3,srt=90,cex=1.2,at=2.5)
      nb1=length(which(data[,1]=="VUS"))
      nb2=length(which(data[,1]=="CON"))
      text(max(maximum,1)/2,3-beforelegend*1.5,paste("                Missense VUS (N=",nb1,") and CI (N=",nb2,")",sep=""),adj=0.5,cex=1.2)
      nb1=length(which(data2[,1]=="VUS"))
      nb2=length(which(data2[,1]=="CON"))
      text(max(maximum,1)/2,2+beforelegend*1.5,paste("LoF VUS (N=",nb1,") and CI (N=",nb2,")",sep=""),adj=0.5,cex=1.2)
      
    }
    
    # plot 1
    {
      PLP=which(data[,1]=="PLP")
      BLB=which(data[,1]=="BLB")
      gnohigh=which(data[,1]=="gnohigh")
      maxaf=10^max(data[PLP,14],-10)
      if(maxaf>0.000001){
        text(max(maximum,1)*legendxfactor*1.02,5.1,paste("Maximal AF in gnomAD = ",format(maxaf*100,digits=2)," %",sep=""),adj=0,cex=1,font=2)
      } else {
        text(max(maximum,1)*legendxfactor*1.02,5.1,"Maximal AF in gnomAD is never seen",adj=0,cex=1,font=2)
      }
      
      if(length(PLP)==0){
        text(max(maximum,1)/2,5.5+0.2,"No PLP missense variant reported",adj=0.5,cex=1.5,col=2,font=3)
      }
      if(length(BLB)==0 & length(gnohigh)==0){
        text(max(maximum,1)/2,5.5-0.2,"No BLB missense variant reported",adj=0.5,cex=1.5,col=2,font=3)
      }
      
      b1=5.5+demigene
      b2=5.5-demigene
      
      if (length(PLP)>0){
        tab=table(data[PLP,18])
        loc=as.numeric(names(tab))
        numb=as.vector(tab)
        for (i in 1:length(loc)){
          lines(c(loc[i],loc[i]),c(b1,b1+step1*numb[i]),col='darkgrey')
        }
        points(loc,b1+step1*numb,col='orange',pch=19)
        #text(loc,b1+step1*numb,loc,adj=0.5,cex=0.2)
      }
      
      lev=sort(unique(c(data[gnohigh,18],data[BLB,18])))
      
      if(length(lev)>0){
        tab=table(factor(data[gnohigh,18],levels=lev))
        loc=as.numeric(names(tab))
        numb=as.vector(tab)
        tab2=table(factor(data[BLB,18],levels=lev))
        numb2=as.vector(tab2)
        
        for (i in 1:length(loc)){
          if(numb[i]>0){
            lines(c(loc[i],loc[i]),c(b2,b2-step1*numb[i]),col='darkgrey')
            points(loc[i],b2-step1*numb[i],col=4,pch=19)
          }
          if(numb2[i]>0){
            lines(c(loc[i],loc[i]),c(b2,b2-step1*numb2[i]),col='darkgrey')
            points(loc[i],b2-step1*numb2[i],col=3,pch=19)
          }
        }
        #text(loc,b2-step1*numb,loc,adj=0.5,cex=0.2)
        #text(loc,b2-step1*numb2,loc,adj=0.5,cex=0.2)
      }
 
    }
    
    
    # plot 2
    {
      PLP=which(data2[,1]=="PLP")
      BLB=which(data2[,1]=="BLB")
      gnohigh=which(data2[,1]=="gnohigh")
      if(length(PLP)>0){
        maxaf=10^max(data2[PLP,14],-10)
      }
      if(length(PLP)==0){
        maxaf=0
      }
      if(maxaf>0.000001){
        text(max(maximum,1)*legendxfactor*1.02,4.1,paste("Maximal AF in gnomAD = ",format(maxaf*100,digits=2)," %",sep=""),adj=0,cex=1,font=2)
      } else {
        text(max(maximum,1)*legendxfactor*1.02,4.1,"Maximal AF in gnomAD is never seen",adj=0,cex=1,font=2)
      }
      
      if(length(PLP)==0){
        text(max(maximum,1)/2,4.5+0.2,"No PLP LoF variant reported",adj=0.5,cex=1.5,col=2,font=3)
      }
      if(length(BLB)==0 & length(gnohigh)==0){
        text(max(maximum,1)/2,4.5-0.2,"No BLB LoF variant reported",adj=0.5,cex=1.5,col=2,font=3)
      }
  
      b1=4.5+demigene
      b2=4.5-demigene

      if (length(PLP)>0){
        tab=table(data2[PLP,18])
        loc=as.numeric(names(tab))
        numb=as.vector(tab)
        for (i in 1:length(loc)){
          lines(c(loc[i],loc[i]),c(b1,b1+step2*numb[i]),col='darkgrey')
        }
        points(loc,b1+step2*numb,col=2,pch=19)
        #text(loc,b1+step2*numb,loc,adj=0.5,cex=0.2)
      }
      
      
      if (length(gnohigh)>0){
        tab=table(data2[gnohigh,18])
        loc=as.numeric(names(tab))
        numb=as.vector(tab)
        for (i in 1:length(loc)){
          lines(c(loc[i],loc[i]),c(b2,b2-step2*numb[i]),col='darkgrey')
        }
        points(loc,b2-step2*numb,col='darkviolet',pch=19)
        #text(loc,b2-step2*numb,loc,adj=0.5,cex=0.2)
      }
      
      if (length(BLB)>0){
        tab=table(data2[BLB,18])
        loc=as.numeric(names(tab))
        numb=as.vector(tab)
        for (i in 1:length(loc)){
          lines(c(loc[i],loc[i]),c(b2,b2-step2*numb[i]),col='darkgrey')
        }
        points(loc,b2-step2*numb,col='darkgreen',pch=19)
        #text(loc,b2-step2*numb,loc,adj=0.5,cex=0.2)
      }
      
    }
    
    
    # # plot 4
    {
      data=databoth[which(databoth[,1]=="VUS" | databoth[,1]=="CON" ),]
      if(dim(data)[1]==0){
        text(max(maximum,1)/2,2.5+0.2,"No missense VUS and conflicting variant reported",adj=0.5,cex=1.5,col=2,font=3)
        text(max(maximum,1)/2,2.5-0.2,"No LoF VUS and conflicting variant reported",adj=0.5,cex=1.5,col=2,font=3)
      }
      if(dim(data)[1]>0){
        b1=2.5+demigene
        b2=2.5-demigene
        VUSm=which(data[,1]=="VUS" & data[,4]=="nonsynonymous SNV")
        CONm=which(data[,1]=="CON" & data[,4]=="nonsynonymous SNV")
        VUSl=which(data[,1]=="VUS" & (data[,4]=="frameshift deletion" | data[,4]=="frameshift insertion" | data[,4]=="stopgain"))
        CONl=which(data[,1]=="CON" & (data[,4]=="frameshift deletion" | data[,4]=="frameshift insertion" | data[,4]=="stopgain"))
        
        if((length(VUSm)+length(CONm))==0){
          text(max(maximum,1)/2,2.5+0.2,"No missense VUS and conflicting variant reported",adj=0.5,cex=1.5,col=2,font=3)
        }
        if((length(VUSl)+length(CONl))==0){
          text(max(maximum,1)/2,2.5-0.2,"No LoF VUS and conflicting variant reported",adj=0.5,cex=1.5,col=2,font=3)
        }
        
        lev=sort(unique(c(data[VUSm,18],data[CONm,18])))
        
        if(length(lev)>0){
          tab=table(factor(data[VUSm,18],levels=lev))
          loc=as.numeric(names(tab))
          numb=as.vector(tab)
          tab2=table(factor(data[CONm,18],levels=lev))
          numb2=as.vector(tab2)
          
          for (i in 1:length(loc)){
            if(numb[i]>0){
              lines(c(loc[i],loc[i]),c(b1,b1+step4*numb[i]),col='darkgrey')
              points(loc[i],b1+step4*numb[i],col=4,pch=19)
            }
            if(numb2[i]>0){
              lines(c(loc[i],loc[i]),c(b1,b1+step4*numb2[i]),col='darkgrey')
              points(loc[i],b1+step4*numb2[i],col=3,pch=19)
            }
          }
          #text(loc,b1+step4*numb,loc,adj=0.5,cex=0.2)
          #text(loc,b1+step4*numb2,loc,adj=0.5,cex=0.2)
        }
        
        lev=sort(unique(c(data[VUSl,18],data[CONl,18])))
        
        if(length(lev)>0){
          tab=table(factor(data[VUSl,18],levels=lev))
          loc=as.numeric(names(tab))
          numb=as.vector(tab)
          tab2=table(factor(data[CONl,18],levels=lev))
          numb2=as.vector(tab2)
          
          for (i in 1:length(loc)){
            if(numb[i]>0){
              lines(c(loc[i],loc[i]),c(b2,b2-step4*numb[i]),col='darkgrey')
              points(loc[i],b2-step4*numb[i],col='orange',pch=19)
            }
            if(numb2[i]>0){
              lines(c(loc[i],loc[i]),c(b2,b2-step4*numb2[i]),col='darkgrey')
              points(loc[i],b2-step4*numb2[i],col=2,pch=19)
            }
          }
          #text(loc,b2-step4*numb,loc,adj=0.5,cex=0.2)
          #text(loc,b2-step4*numb2,loc,adj=0.5,cex=0.2)
        }
        
      }
    }
    
    
    # plot 3
    {
      if(dim(cons)[1]>50){
        posi=1:length(consw)
        rem=which(posi>=zoom[1] & posi<=zoom[2])
        posi=(posi-zoom[1])/(zoom[2]-zoom[1])*maximum
        lines(posi[rem],consw[rem],col=3,lwd=3)
      }
      
      # PROTEIN DOMAINS
      
      # take id of protein
      {
        defaultW <- getOption("warn") 
        options(warn = -1) 
        id1=strsplit(as.character(isoform), split="\\.")[[1]][1]
        options(warn = defaultW)
        id2=refseq[which(refseq[,3]==id1),1]
        if (length(id2)>1){
          id3=""
          for (iso in 1:length(id2)){
            id3=c(id3,as.character(uniprot[which(as.character(uniprot[,1])==id2[iso]),2]))
          }
          if(length(which(id3==""))>0){
            id3=id3[-which(id3=="")]
            id3=id3[1]
          }
        } else {
          id3=as.character(uniprot[which(as.character(uniprot[,1])==id2),2])
        }
      }
      
      #id3=""
      
      # plot domains and region
      if(length(id3)>0){if(id3!="" & is.na(id3)=="FALSE"){
        
        # download data about the protein
        print(id3)
        rel_json <- drawProteins::get_features(id3)
        prot <- drawProteins::feature_to_dataframe(rel_json)
        prot=prot[sort(prot$begin, index.return=TRUE)$ix,]
        
        # computing size of annotations
        l=1:length(prot$type)
        for (g in 1:length(l)){
          l[g]=prot$end[g]-prot$begin[g]
        }
        
        # for domains and regions
        
        all=which(l>0 & (prot$type=="REGION" | prot$type=="DOMAIN" | prot$type=="ZN_FING" | prot$type=="REPEAT" | prot$type=="CHAIN" | prot$type=="TRANSIT" | prot$type=="DNA_BIND" | prot$type=="TRANSMEM" | prot$type=="TOPO_DOM") )
        prot[all,2]=gsub(" \\d+","",prot[all,2])
        prot[,2]=sapply(strsplit(as.character(prot[,2]), ";"), `[`, 1)
        if(length(unique(prot[all,2]))>=21){
         all=all[1:20] 
        }
        coldom=GetColors(n = length(unique(prot[all,2])), scheme = "discrete rainbow")
        
        
        
        # CHAINS
        dom1=which(l>0 &  prot$type=="CHAIN")
        if(length(dom1)>0){
          uni=unique(prot$description[dom1])
          ndom1=length(unique(prot$description[dom1]))
          cols=coldom[1:ndom1]
          for (d in 1:length(uni)){
            a=which(prot$description==uni[d])
            for (e in 1:length(a)){
              prot$description[a[e]]="Mature chain"
              
              pos1=min(prot$begin[a[e]],max(maximum,1))
              pos2=min(prot$end[a[e]],max(maximum,1))
              if(pos1<zoom[1]){pos1=zoom[1]}
              if(pos2>zoom[2]){pos2=zoom[2]}
              
              if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
                pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
                pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
                rect(pos1,    3.5-demigene,   pos2,   3.5+demigene,   col=cols[d])
              }
            }
          }
        } else {
          ndom1=0
        }
        
          # domains
        dom2=which(l>0 & (prot$type=="DOMAIN" | prot$type=="ZN_FING" | prot$type=="DNA_BIND"))
        prot[which(prot$type=="ZN_FING"),2]=paste("ZN_FING ",prot[which(prot$type=="ZN_FING"),2],sep="")
        if(length(dom2)>0){
          uni=unique(prot$description[dom2])
          ndom2=length(unique(prot$description[dom2]))
          cols=coldom[(ndom1+1):(ndom1+ndom2)]
          for (d in 1:length(uni)){
            a=which(prot$description==uni[d])
            for (e in 1:length(a)){
              pos1=min(max(maximum,1),prot$begin[a[e]],na.rm=T)
              pos2=min(prot$end[a[e]],max(maximum,1),na.rm=T)
              if(pos1<zoom[1]){pos1=zoom[1]}
              if(pos2>zoom[2]){pos2=zoom[2]}
              
              if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
                pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
                pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
                rect(pos1,   3.5-demigene+0.015,  pos2 ,   3.5+demigene-0.015,col=cols[d])
              }
            }
          }
          
        } else {
          ndom2=0
        }
        
        ### FOR REGIONs
        dom3=which(l>0 & (prot$type=="REGION"| prot$type=="REPEAT" | prot$type=="TRANSIT" | prot$type=="TRANSMEM" | prot$type=="TOPO_DOM"))
        if(length(dom3)>0){
          uni=unique(prot$description[dom3])
          ndom3=length(unique(prot$description[dom3]))
          cols=coldom[(ndom1+ndom2+1):(ndom1+ndom2+ndom3)]
          for (d in 1:length(uni)){
            a=which(prot$description==uni[d])
            for (e in 1:length(a)){
              
              pos1=min(prot$begin[a[e]],max(maximum,1))
              pos2=min(prot$end[a[e]],max(maximum,1))
              if(pos1<zoom[1]){pos1=zoom[1]}
              if(pos2>zoom[2]){pos2=zoom[2]}
              
              if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
                pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
                pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
                rect(pos1,   3.5-demigene+0.03,  pos2  ,   3.5+demigene-0.03,col=cols[d])
              }
            }
          }
        }
        if(length(all)>0){
          text(0,3.35,"Mature chains (large), domains (medium) and regions (narrow) from Uniprot (>=10aa)",adj=0,cex=1.2)
          legend(0,3.28,unique(prot$description[c(dom1,dom2,dom3)]),fill=coldom,col=coldom,cex=0.5,ncol=3,bg='white')
        }
      }}
    }
    

    # clusters
    {
    clustPLP=0
    clustBLB=0
    
    clup=regionsp[which(regionsp[,2]==isoform),]
    if(length(which(regionsp[,2]==isoform))>1){
      clustPLP=1
      for (i in 1:dim(clup)[1]){
        pos1=as.numeric(clup[i,3])
        pos2=as.numeric(clup[i,4])
        if(pos1<zoom[1]){pos1=zoom[1]}
        if(pos2>zoom[2]){pos2=zoom[2]}
        if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
          pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
          pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
          rect(pos1,5.5+demigene-0.015,pos2,5.5+demigene+0.015,col='orange')
        }
        if(pos1==zoom[1] & pos2==zoom[2]){
          rect(0,5.5+demigene-0.015,maximum,5.5+demigene+0.015,col='orange')
        }
      }
    }
    if(length(which(regionsp[,2]==isoform))==1){
      clustPLP=1
      pos1=as.numeric(clup[3])
      pos2=as.numeric(clup[4])
      if(pos1<zoom[1]){pos1=zoom[1]}
      if(pos2>zoom[2]){pos2=zoom[2]}
      if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
        pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
        pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
        rect(pos1,5.5+demigene-0.015,pos2,5.5+demigene+0.015,col='orange')
      }
      if(pos1==zoom[1] & pos2==zoom[2]){
        rect(0,5.5+demigene-0.015,maximum,5.5+demigene+0.015,col='orange')
      }
    }
    
    club=regionsb[which(regionsb[,2]==isoform),]
    if(length(which(regionsb[,2]==isoform))>1){
      clustBLB=1
      for (i in 1:dim(club)[1]){
        pos1=as.numeric(club[i,3])
        pos2=as.numeric(club[i,4])
        if(pos1<zoom[1]){pos1=zoom[1]}
        if(pos2>zoom[2]){pos2=zoom[2]}
        if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
          pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
          pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
          rect(pos1,5.5-demigene-0.015,pos2,5.5-demigene+0.015,col='turquoise3')
        }
        if(pos1==zoom[1] & pos2==zoom[2]){
          rect(0,5.5-demigene-0.015,maximum,5.5-demigene+0.015,col='turquoise3')
        }
      }
    }
    if(length(which(regionsb[,2]==isoform))==1){
      clustBLB=1
      pos1=as.numeric(club[3])
      pos2=as.numeric(club[4])
      if(pos1<zoom[1]){pos1=zoom[1]}
      if(pos2>zoom[2]){pos2=zoom[2]}
      if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
        pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
        pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
        rect(pos1,5.5-demigene-0.015,pos2,5.5-demigene+0.015,col='turquoise3')
      }
      if(pos1==zoom[1] & pos2==zoom[2]){
        rect(0,5.5-demigene-0.015,maximum,5.5-demigene+0.015,col='turquoise3')
      }
    }
    
    # legend of clusters
    
    demiclust=0.015
    end=maximum/35
    beg=maximum/29
    
    if(clustBLB==1 & clustPLP==1){
      rect(0,5.05-demiclust,end,5.05+demiclust,col='turquoise3')
      rect(0,5.11-demiclust,end,5.11+demiclust,col='orange')
      text(beg,5.11,"Significant PLP clustering",adj=0,cex=0.9)
      text(beg,5.05,"Significant BLB clustering",adj=0,cex=0.9)
    }
    if(clustBLB==1 & clustPLP==0){
      rect(0,5.075-demiclust,end,5.075+demiclust,col='turquoise3')
      text(beg,5.075,"Significant BLB clustering",adj=0,cex=0.9)
    }
    if(clustBLB==0 & clustPLP==1){
      rect(0,5.075-demiclust,end,5.075+demiclust,col='orange')
      text(beg,5.075,"Significant PLP clustering",adj=0,cex=0.9)
    }
    if(clustBLB==0 & clustPLP==0){
      text(0,5.07,"No significant clustering detected",adj=0,cex=0.9)
    }
  
    }
    
    }
    

    time2=Sys.time()
    print(time2-time1)
    
  })
  
  output$plot1 <- renderPlot({
    print("Printing plot")
    print(plotInput())
  })
  
  plotInput2 <- reactive({
    print("Doing plot2")
    time1=Sys.time()

    dataallgenes=variables$dataallgenes2
    bla=dataInput()
    isoform=bla[1]
    data=dataallgenes[which(dataallgenes[,2]==isoform),]

    if(dim(data)[1]>0){
      
      mut1=variables$mut2
      
      bla=dataInput()
      
      isoform=bla[1]
      zoom=c(as.numeric(bla[2]),as.numeric(bla[3]))
      constrack=bla[4]
      conswindow=as.numeric(bla[5])
      cchange=bla[6]
      pchange=bla[7]
      radio=bla[8]
      
      if(is.na(pchange)){pchange="-"}
      if(is.na(cchange)){cchange="-"}
      if(is.na(radio)){radio="DNA"}
      
      if(pchange!="-" & radio=="Protein"){
        cchange=mut1[which(mut1[,1]==isoform & mut1[,4]==pchange)[1],2]
      }
      
      if(pchange=="-" & radio=="Protein"){
        cchange="-"
      }
      
      if(cchange!="-" & radio=="DNA"){
        pchange=mut1[which(mut1[,1]==isoform & mut1[,2]==cchange),4]
      }
      
      if(cchange=="-" & radio=="DNA"){
        pchange="-"
      }
  
      load(file="PLP-BLB.RData")
      
      par(mfrow=c(3,1))
      hist(PLP,n=50,xlab="MutScore",cex.main=1.7,cex.lab=1.5,cex.axis=1.3,col='red',main="MutScore distribution for PLP missense variants from ClinVar",yaxt='n',ylab="")
      title(ylab="Frequency", line=(2), cex.lab=1.5)
      axis(2,at=c(0,3000,6000,9000,12000),labels=c(0,3000,6000,9000,12000),las=2,cex.axis=1.3,pos=0)
      
      if(cchange!="-" & pchange!="-"){
        muto=mut1[which(mut1[,2]==cchange & mut1[,1]==isoform),6]
        abline(v=muto,lwd=2,col='green')
        legend(0.5, 13000, legend=paste(dataallgenes[which(dataallgenes[,2]==isoform)[1],5],", ",isoform,":",cchange,", ",pchange,"    MutScore = ",muto,sep=""),col=c("green"), lty=1, cex=1.5, xjust=0.5,lwd=2,bg="white")
        
        n=50
        a=hist(PLP,n=n,plot=F)
        b=hist(BLB,n=n,plot=F)
        m=which.min(abs(a$mids-muto))[1]
        proba=round(a$density[m]/(a$density[m]+b$density[m])*100,digits=1)
        text(0.5,6000,paste("",proba,"% of variants with similar scores are PLPs",sep=""),cex=1.3)
        proba2=round(length(which(PLP>muto))/length(PLP)*100,digits=1)
        text(0.5,3500,paste("",proba2,"% of PLPs have a score higher than ",muto,sep=""),cex=1.3)
      }
      
      hist(BLB,n=50,xlab="MutScore",cex.main=1.7,cex.lab=1.5,cex.axis=1.3,col='blue',main="MutScore distribution for BLB missense variants from ClinVar",yaxt='n',ylab="")
      title(ylab="Frequency", line=(2), cex.lab=1.5)
      axis(2,at=c(0,1000,2000,3000),labels=c(0,1000,2000,3000),las=2,cex.axis=1.3,pos=0)
      if(cchange!="-" & pchange!="-"){
        abline(v=muto,lwd=2,col='green')
        legend(0.5, 3900, legend=paste(dataallgenes[which(dataallgenes[,2]==isoform)[1],5],", ",isoform,":",cchange,", ",pchange,"    MutScore = ",muto,sep=""),col=c("green"), lty=1, cex=1.5, xjust=0.5,lwd=2,bg="white")
        text(0.5,1800,paste("",100-proba,"% of variants with similar scores are BLBs",sep=""),cex=1.3)
        proba2=round(length(which(BLB<muto))/length(BLB)*100,digits=1)
        text(0.5,1100,paste("",proba2,"% of BLBs have a score lower than ",muto,sep=""),cex=1.3)
      }

    }
    
    time2=Sys.time()
    print(time2-time1)
    
  
  })
  
  output$plot2 <- renderPlot({
    print("Printing plot2")
    print(plotInput2())
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { 
      
      isoform=input$isoform
      
      cchange=input$cpos
      
      pchange=input$ppos
      
      radio=input$radio

      if(pchange!="-" & radio=="Protein"){
        cchange=mut1[which(mut1[,1]==isoform & mut1[,4]==pchange)[1],2]
      }
      
      if(pchange=="-" & radio=="Protein"){
        cchange="-"
      }
      
      if(cchange!="-" & radio=="DNA"){
        pchange=mut1[which(mut1[,1]==isoform & mut1[,2]==cchange),4]
      }
      
      if(cchange=="-" & radio=="DNA"){
        pchange="-"
      }
      
      if(input$cpos!="-") {return(paste("MutLand_",input$gene,"_",input$isoform,"_",cchange,"_",pchange,'.pdf', sep=''))}
      if(input$cpos=="-") {return(paste("MutLand_",input$gene,"_",input$isoform,'.pdf', sep=''))}
      },
    content = function(file){
      pdf(file, height=12,width=17)
      
      print("DOING PLOT for PDF")
      
      time1=Sys.time()
      mut1=variables$mut2
      data1=variables$data12
      dataallgenes=variables$dataallgenes2
      scores=variables$scores2
      
      bla=dataInput()
      
      isoform=bla[1]
      zoom=c(as.numeric(bla[2]),as.numeric(bla[3]))
      constrack=bla[4]
      conswindow=as.numeric(bla[5])
      cchange=bla[6]
      pchange=bla[7]
      radio=bla[8]
      
      data=dataallgenes[which(dataallgenes[,2]==isoform),]
      
      print(paste("Check plot: ",pchange,", ",cchange,", ",radio,", ",dim(data)[1],sep=""))
      print(isoform)
      
      if(is.na(pchange)==FALSE & is.na(cchange)==FALSE & is.na(radio)==FALSE & dim(data)[1]>0 & zoom[2]>1){
        
        if(is.na(pchange)){pchange="-"}
        if(is.na(cchange)){cchange="-"}
        if(is.na(radio)){radio="DNA"}
        
        if(pchange!="-" & radio=="Protein"){
          cchange=mut1[which(mut1[,1]==isoform & mut1[,4]==pchange)[1],2]
        }
        
        if(pchange=="-" & radio=="Protein"){
          cchange="-"
        }
        
        if(cchange!="-" & radio=="DNA"){
          pchange=mut1[which(mut1[,1]==isoform & mut1[,2]==cchange),4]
        }
        
        if(cchange=="-" & radio=="DNA"){
          pchange="-"
        }
        
        maximum=as.numeric(maxiso[which(maxiso[,1]==isoform),2])
        
        data=dataallgenes[which(dataallgenes[,2]==isoform),]
        
        gene=as.character(data[1,5])
        
        # other types of variants
        datanc1=dataallgenes[which(dataallgenes[,5]==gene),]
        datanc=datanc1[which(datanc1[,4]=="intronic" | datanc1[,4]=="UTR3" | datanc1[,4]=="UTR5" | datanc1[,4]=="splicing"),]
        data=data[which(data[,4]=="frameshift deletion" | data[,4]=="frameshift insertion" | data[,4]=="stopgain" | data[,4]=="nonsynonymous SNV" | data[,4]=="synonymous SNV"),]
        
        data=data[which(data$pos>=zoom[1] & data$pos<=zoom[2]),]
        data$pos=round((data$pos-zoom[1])/(zoom[2]-zoom[1])*maximum)
        
        # choosing frequent gnomAD variants for ClinVar plot
        data[which(data[,14]==0),14]=(-6)
        data[,14]=as.numeric(data[,14])
        data[which(data[,15]=="-7"),1]="gnohigh"
        
        # separating missenses and LoF
        dataall=data
        data2=data[which(data[,4]=="frameshift deletion" | data[,4]=="frameshift insertion" | data[,4]=="stopgain"),]
        data=data[which(data[,4]=="nonsynonymous SNV"),]
        databoth=rbind(data,data2)
        
        # computing steps for graphics of clinvar and gnomad
        {
          # for plot 1
          
          PLP=which(data[,1]=="PLP")
          BLB=which(data[,1]=="BLB")
          gnohigh=which(data[,1]=="gnohigh")
          m=max(table(data[PLP,18]),table(data[BLB,18]),table(data[gnohigh,18]))
          step1=demiinterval/m
          m1=m
          if(m1>100000){m1=1}
          if(m1<1){m1=1}
          
          # for plot 2
          PLP=which(data2[,1]=="PLP")
          BLB=which(data2[,1]=="BLB")
          gnohigh=which(data2[,1]=="gnohigh")
          m=0
          if (dim(data2)[1]>0){
            m=max(table(data2[PLP,18]),table(data2[BLB,18]),table(data2[gnohigh,18]))
          } else {
            m=1
          }
          step2=demiinterval/m
          m2=m
          step2=min(step1,step2)
          m2=max(m1,m2)
          if(m2>100000){m2=1}
          if(m2<1){m2=1}
          
          # for plot 3
          #conservation
          nb=3
          if(constrack=="GERP++RS"){nb=3}
          if(constrack=="PhyloP100way"){nb=4}
          if(constrack=="PhastCons100way"){nb=5}
          if(constrack=="pext"){nb=6}
          if(constrack=="MutScore"){nb=8}
          
          size=max(maximum,1)
          cons=scores[which(scores[,1]==isoform),]
          
          a0=as.matrix(cons[,c(7,..nb)])
          a1=1:size
          a1[1:size]=NA
          for (i in 1:size){
            a1[i]=mean(a0[which(a0[,1]==i),2])
          }
          
          consw=runmean(a1,conswindow)
          #window(as.matrix(cons[,c(7,..nb)]),start=1,end=size,deltat=conswindow)
          
          ma=ceiling(max(consw,na.rm=T))
          mi=floor(min(0,min(consw,na.rm=T)))
          if(mi==ma){mi=mi-1}
          m3=ma-mi
          step3=demiinterval/m3
          fact=1/step3
          consw=3.5+demigene+(consw-mi)/fact
          
          # for plot 4
          VUS=which(data[,1]=="VUS")
          CON=which(data[,1]=="CON")
          m=0
          if (length(VUS)+length(CON)>0){
            m=max(table(data[VUS,18]),table(data[CON,18]))
          } else {
            m=1
          }
          mtemp=m
          m=0
          VUS=which(data2[,1]=="VUS")
          CON=which(data2[,1]=="CON")
          
          if (dim(data2)[1]>0){
            m=max(table(data2[VUS,18]),table(data2[CON,18]))
          } else {
            m=1
          }
          m4=max(m,mtemp)
          step4=demiinterval/m4
          
          if(m3>100000){m3=1}
          if(m3<1){m3=1}
          
          if(m4>100000){m4=1}
          if(m4<1){m4=1}
          
        }
        print(m1)
        print(m2)
        print(m3)
        print(m4)
        
        defaultW <- getOption("warn") 
        options(warn = -1) 
        nm=strsplit(as.character(isoform),"\\.")[[1]][1]
        options(warn = defaultW)
        enst=conv[which(conv[,4]==nm),2][1]
        
        #layout(matrix(c(1,1,1,1,2,3), nrow = 6, ncol = 1, byrow = TRUE))
        # begin to plot
        plot(0,1000,ylim=c(2,6),col=2,xlim=c(0,max(as.numeric(maxiso[which(maxiso[,1]==isoform),2]),1)*1.32), xlab="",ylab="",
             main=bquote(bold("MutLand plot for " ~ bolditalic(.(gene)) ~ " (" ~ .(as.character(isoform)) ~ " - " ~ .(as.character(enst)) ~ ")")),yaxt='n',xaxt='n',cex.main=2,yaxs="i",cex.lab=2.5,font.main=2)
        mtext("Amino acid position",side=1,line=3,cex=1.5,at=max(maximum,1)/2,adj=0.5)
        
        if(cchange!="-"){
          position=mut1[which(mut1[,1]==isoform & mut1[,2]==cchange),5]
          posiline=round((position-zoom[1])/(zoom[2]-zoom[1])*maximum)
          if (position>0){abline(v=posiline,col='red')}
        }
        
        #Grids and rectangle
        {
          # vertical grid
          # choosing size between x-axis labels
          po=c(10,50,100,200,500,1000,2000,4000,6000)
          a=po[which.min(abs(max(zoom[2]-zoom[1],1)/po-10))]
          b=0:20
          ax2=a*b+round(zoom[1]/a)*a
          
          ax=(ax2-zoom[1])/(zoom[2]-zoom[1])*maximum
          
          rem=which(ax2<=zoom[2] & ax2>=zoom[1])
          ax2=ax2[rem]
          ax=ax[rem]
          
          print("Debug")
          print(zoom[1])
          print(zoom[2])
          print(maximum)
          
          axis(1,at=ax,labels=ax2,cex.axis=1.5)
          for (i in 1:length(ax)){
            abline(v=ax[i],lty=3,col='lightgray')
          }
          
          # horizontal grids
          y=c()
          lab=c()
          for (i in 1:m1){
            y=c(y,5.5-demigene-step1*i,5.5+demigene+step1*i)
            lab=c(lab,i,i)
          }
          for (i in 1:m2){
            y=c(y,4.5-demigene-step2*i,4.5+demigene+step2*i)
            lab=c(lab,i,i)
          }
          for (i in 0:m3){
            y=c(y,3.5+demigene+step3*i)
            lab=c(lab,i+mi)
          }
          for (i in 1:m4){
            y=c(y,2.5-demigene-step4*i,2.5+demigene+step4*i)
            lab=c(lab,i,i)
          }
          axis(2, at = y,lab,cex.axis=0.6)
          for (i in 1:length(y)){
            abline(h=y[i],lty=3,col='lightgray')
          }
          abline(h=1)
          abline(h=2)
          abline(h=3)
          abline(h=4)
          abline(h=5)
          
          # gene rectangles
          rect(0,5.5-demigene,max(maximum,1),5.5+demigene,col='lightgray')
          rect(0,4.5-demigene,max(maximum,1),4.5+demigene,col='lightgray')
          rect(0,3.5-0.03,max(maximum,1),3.5+0.03,col='lightgray')
          rect(0,2.5-demigene,max(maximum,1),2.5+demigene,col='lightgray')
          rect(0,1.5-demigene,max(maximum,1),1.5+demigene,col='lightgray')
          rect(0,0.5-demigene,max(maximum,1),0.5+demigene,col='lightgray')
          
        }
        
        # legends
        {
          # plot 1
          text(max(maximum,1)*legendxfactor*(-0.03),6-beforelegend*1.5,"ClinVar missense variants",adj=0,cex=1.5,font=2)
          legend(max(maximum,1)*legendxfactor,6-beforelegend*1,
                 c("ClinVar* missense PLP variants","ClinVar* missense BLB variants","Frequent** gnomAD missense variants "),
                 col=c('orange',3,4),pch=c(19,19,19),lty=NA,lwd=4,cex=1,bg='white', pt.cex=1.2)
          text(max(maximum,1)*legendxfactor*1.02,2.1,"** with AF > max(AF of PLPs)]",adj=0,cex=1.2)
          text(max(maximum,1)*legendxfactor*1.02,2.2, "* ClinVar version of 21.11.2020",adj=0,cex=1.2)
          mtext("Missense variants per aa",side=2,line=3,srt=90,cex=1.2,at=5.5)
          nb1=length(which(data[,1]=="PLP"))
          text(max(maximum,1)/2,6-beforelegend*1.5,paste("PLP missense variants from ClinVar (N=",nb1,")",sep=""),adj=0.5,cex=1.2)
          nb1=length(which(data[,1]=="BLB"))
          text(max(maximum,1)/2,5+beforelegend*1.5,paste("BLB missense variants from ClinVar (N=",nb1,")",sep=""),adj=0.5,cex=1.2)
          
          nb1=length(which(datanc[,4]=="intronic" & datanc[,1]=="PLP"))
          nb2=length(which(datanc[,4]=="UTR5" & datanc[,1]=="PLP"))
          nb3=length(which(datanc[,4]=="UTR3" & datanc[,1]=="PLP"))
          nb4=length(which(dataall[,4]=="synonymous SNV" & dataall[,1]=="PLP"))
          if(nb1+nb2+nb3+nb4>0){
            text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*9,"In ClinVar there are also:",adj=0,cex=1.2)
            bot=11.5
            if(nb1>0){
              text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*bot,paste("- ",nb1, " intronic PLP variant(s)",sep=""),adj=0,cex=1)
              bot=bot+2
            }
            if(nb2>0){
              text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*bot,paste("- ",nb2, " 5'UTR PLP variant(s)",sep=""),adj=0,cex=1)
              bot=bot+2
            }
            if(nb3>0){
              text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*bot,paste("- ",nb3, " 3'UTR PLP variant(s)",sep=""),adj=0,cex=1)
              bot=bot+2
            }
            if(nb4>0){
              text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*bot,paste("- ",nb4, " synonymous PLP variant(s)",sep=""),adj=0,cex=1)
              bot=bot+2
            }
          } else {
            text(max(maximum,1)*legendxfactor*1.02,6-beforelegend*9,"No non-coding or synonymous PLP",adj=0,cex=1)
          }
          
          # exons
          
          exons=unique(cons[,c(2,7)])
          nbexons=max(exons[,1])
          minexon=min(exons[,1])
          posi3=1
          text(0,5.5,"Exons ",adj=1,cex=0.8)
          for (i in minexon:nbexons){
            
            if(length(which(exons[,1]==i))>0){
              
              posi=max(exons[which(exons[,1]==i),2])
              
              posi2=round((posi-zoom[1])/(zoom[2]-zoom[1])*maximum)
              if(posi2>=0 & posi2<=maximum){
                rect(posi2,5.5-demigene,posi2,5.5+demigene)
              }
              if((posi2-posi3)>max(maximum,1)/50 & (posi2-posi3)/2+posi3>maximum*0.01 & (posi2-posi3)/2+posi3<maximum*0.99){
                text((posi2-posi3)/2+posi3,5.5,i,adj=0.5,cex=0.8)
              }
              if(posi2=="-Inf"){posi2=1}
              posi3=posi2
              
            }
          }
          
          
          # plot 2
          text(max(maximum,1)*legendxfactor*(-0.03),5-beforelegend*1.5,"Clinvar LoF variants",adj=0,cex=1.5,font=2)
          legend(max(maximum,1)*legendxfactor,5-beforelegend*1,
                 c("ClinVar* LoF variants","ClinVar* LoF BLB variants","Frequent** gnomAD LoF variants"),
                 col=c(2,'darkgreen','darkviolet'),pch=c(19,19,19),lty=NA,lwd=4,cex=1,bg='white',pt.cex=1.2)
          mtext("LoF variants per aa",side=2,line=3,srt=90,cex=1.2,at=4.5)
          nb2=length(which(data2[,1]=="PLP"))
          text(max(maximum,1)/2,5-beforelegend*1.5,paste("PLP LoF variants from ClinVar (N=",nb2,")",sep=""),adj=0.5,cex=1.2)
          nb2=length(which(data2[,1]=="BLB"))
          text(max(maximum,1)/2,4+beforelegend*1.5,paste("BLB LoF variants from ClinVar (N=",nb2,")",sep=""),adj=0.5,cex=1.2)
          
          nb1=length(which(datanc[,4]=="splicing" & datanc[,1]=="PLP"))
          if (nb1>0){
            text(max(maximum,1)*legendxfactor*1.02,5-beforelegend*9,"In ClinVar there are also:",adj=0,cex=1.2)
            text(max(maximum,1)*legendxfactor*1.02,5-beforelegend*11.5,paste("- ",nb1, " splicing PLP variant(s)",sep=""),adj=0,cex=1)
          } else {
            text(max(maximum,1)*legendxfactor*1.02,5-beforelegend*9,"No splicing PLP reported",adj=0,cex=1)
          }
          
          # exons
          
          exons=unique(cons[,c(2,7)])
          nbexons=max(exons[,1])
          minexon=min(exons[,1])
          posi3=1
          text(0,4.5,"Exons ",adj=1,cex=0.8)
          for (i in minexon:nbexons){
            if(length(which(exons[,1]==i))>0){
              posi=max(exons[which(exons[,1]==i),2])
              posi2=round((posi-zoom[1])/(zoom[2]-zoom[1])*maximum)
              if(posi2>=0 & posi2<=maximum){
                rect(posi2,4.5-demigene,posi2,4.5+demigene)
              }
              if((posi2-posi3)>max(maximum,1)/50 & (posi2-posi3)/2+posi3>maximum*0.01 & (posi2-posi3)/2+posi3<maximum*0.99){
                text((posi2-posi3)/2+posi3,4.5,i,adj=0.5,cex=0.8)
              }
              if(posi2=="-Inf"){posi2=1}
              posi3=posi2
            }
          }
          
          # plot 3
          legend(max(maximum,1)*legendxfactor,4-beforelegend*2,paste(constrack," ","***",sep=""),col=3,lty=1,lwd=4,cex=1,bg='white')
          mtext("Score",side=2,line=3,srt=90,cex=1.2,at=3.65)
          text(max(maximum,1)*legendxfactor*(-0.03),4-beforelegend*1.5,"Protein domains and scores",cex=1.5,font=2,adj=0)
          text(max(maximum,1)*legendxfactor*1.02,4-beforelegend*6.5, paste("*** averaged with a window of size ",conswindow,sep=""),adj=0,cex=0.8)
          
          if(length(which(data[,1]=="PLP"))<10){
            text(max(maximum,1)*legendxfactor*1.02,4-beforelegend*8,"Note: MutScore computed without positional score",adj=0,cex=0.8)
          }
          
          if(cchange!="-"){
            text(max(maximum,1)*legendxfactor,4-beforelegend*13, paste("",cchange," - ",pchange,sep=""),adj=0,cex=1.7,col=2,font=2)
            muto=mut1[which(mut1[,2]==cchange & mut1[,1]==isoform),6]
            text(max(maximum,1)*legendxfactor,4-beforelegend*16, paste("MutScore = ",muto,sep=""),adj=0,cex=1.7,col=2,font=2)
            rect(max(maximum,1)*legendxfactor*0.99, 4-beforelegend*18, max(maximum,1)*legendxfactor*1.3, 4-beforelegend*11,border=2,lwd=2)
          }
          
          # exons
          
          exons=unique(cons[,c(2,7)])
          nbexons=max(exons[,1])
          minexon=min(exons[,1])
          posi3=1
          text(0,2.5,"Exons ",adj=1,cex=0.8)
          for (i in minexon:nbexons){
            if(length(which(exons[,1]==i))>0){
              posi=max(exons[which(exons[,1]==i),2])
              posi2=round((posi-zoom[1])/(zoom[2]-zoom[1])*maximum)
              if(posi2>=0 & posi2<=maximum){
                rect(posi2,2.5-demigene,posi2,2.5+demigene)
              }
              if((posi2-posi3)>max(maximum,1)/50 & (posi2-posi3)/2+posi3>maximum*0.01 & (posi2-posi3)/2+posi3<maximum*0.99){
                text((posi2-posi3)/2+posi3,2.5,i,adj=0.5,cex=0.8)
              }
              if(posi2=="-Inf"){posi2=1}
              posi3=posi2
            }
          }
          
          # plot 6
          text(max(maximum,1)*legendxfactor*(-0.03),3-beforelegend*1.5,"VUS and CI variants",adj=0,cex=1.5,font=2)
          label=c("ClinVar* VUS missense","ClinVar* CI missense","ClinVar* VUS LoF","ClinVar* CI LoF")
          text(max(maximum,1)*legendxfactor*1.02,3-beforelegend*8.2, "* ClinVar version of 21.11.2020",adj=0,cex=0.8)
          legend(max(maximum,1)*legendxfactor,3-beforelegend*4,label,col=c(4,3,'orange',2),lty=NA,lwd=4,pch=c(19,19,19,19),cex=1,bg='white',pt.cex=1.2)
          mtext("Variants per aa",side=2,line=3,srt=90,cex=1.2,at=2.5)
          nb1=length(which(data[,1]=="VUS"))
          nb2=length(which(data[,1]=="CON"))
          text(max(maximum,1)/2,3-beforelegend*1.5,paste("                Missense VUS (N=",nb1,") and CI (N=",nb2,")",sep=""),adj=0.5,cex=1.2)
          nb1=length(which(data2[,1]=="VUS"))
          nb2=length(which(data2[,1]=="CON"))
          text(max(maximum,1)/2,2+beforelegend*1.5,paste("LoF VUS (N=",nb1,") and CI (N=",nb2,")",sep=""),adj=0.5,cex=1.2)
          
        }
        
        # plot 1
        {
          PLP=which(data[,1]=="PLP")
          BLB=which(data[,1]=="BLB")
          gnohigh=which(data[,1]=="gnohigh")
          maxaf=10^max(data[PLP,14],-10)
          if(maxaf>0.000001){
            text(max(maximum,1)*legendxfactor*1.02,5.1,paste("Maximal AF in gnomAD = ",format(maxaf*100,digits=2)," %",sep=""),adj=0,cex=1,font=2)
          } else {
            text(max(maximum,1)*legendxfactor*1.02,5.1,"Maximal AF in gnomAD is never seen",adj=0,cex=1,font=2)
          }
          
          if(length(PLP)==0){
            text(max(maximum,1)/2,5.5+0.2,"No PLP missense variant reported",adj=0.5,cex=1.5,col=2,font=3)
          }
          if(length(BLB)==0 & length(gnohigh)==0){
            text(max(maximum,1)/2,5.5-0.2,"No BLB missense variant reported",adj=0.5,cex=1.5,col=2,font=3)
          }
          
          b1=5.5+demigene
          b2=5.5-demigene
          
          if (length(PLP)>0){
            tab=table(data[PLP,18])
            loc=as.numeric(names(tab))
            numb=as.vector(tab)
            for (i in 1:length(loc)){
              lines(c(loc[i],loc[i]),c(b1,b1+step1*numb[i]),col='darkgrey')
            }
            points(loc,b1+step1*numb,col='orange',pch=19)
            #text(loc,b1+step1*numb,loc,adj=0.5,cex=0.2)
          }
          
          lev=sort(unique(c(data[gnohigh,18],data[BLB,18])))
          
          if(length(lev)>0){
            tab=table(factor(data[gnohigh,18],levels=lev))
            loc=as.numeric(names(tab))
            numb=as.vector(tab)
            tab2=table(factor(data[BLB,18],levels=lev))
            numb2=as.vector(tab2)
            
            for (i in 1:length(loc)){
              if(numb[i]>0){
                lines(c(loc[i],loc[i]),c(b2,b2-step1*numb[i]),col='darkgrey')
                points(loc[i],b2-step1*numb[i],col=4,pch=19)
              }
              if(numb2[i]>0){
                lines(c(loc[i],loc[i]),c(b2,b2-step1*numb2[i]),col='darkgrey')
                points(loc[i],b2-step1*numb2[i],col=3,pch=19)
              }
            }
            #text(loc,b2-step1*numb,loc,adj=0.5,cex=0.2)
            #text(loc,b2-step1*numb2,loc,adj=0.5,cex=0.2)
          }
          
        }
        
        
        # plot 2
        {
          PLP=which(data2[,1]=="PLP")
          BLB=which(data2[,1]=="BLB")
          gnohigh=which(data2[,1]=="gnohigh")
          if(length(PLP)>0){
            maxaf=10^max(data2[PLP,14],-10)
          }
          if(length(PLP)==0){
            maxaf=0
          }
          if(maxaf>0.000001){
            text(max(maximum,1)*legendxfactor*1.02,4.1,paste("Maximal AF in gnomAD = ",format(maxaf*100,digits=2)," %",sep=""),adj=0,cex=1,font=2)
          } else {
            text(max(maximum,1)*legendxfactor*1.02,4.1,"Maximal AF in gnomAD is never seen",adj=0,cex=1,font=2)
          }
          
          if(length(PLP)==0){
            text(max(maximum,1)/2,4.5+0.2,"No PLP LoF variant reported",adj=0.5,cex=1.5,col=2,font=3)
          }
          if(length(BLB)==0 & length(gnohigh)==0){
            text(max(maximum,1)/2,4.5-0.2,"No BLB LoF variant reported",adj=0.5,cex=1.5,col=2,font=3)
          }
          
          b1=4.5+demigene
          b2=4.5-demigene
          
          if (length(PLP)>0){
            tab=table(data2[PLP,18])
            loc=as.numeric(names(tab))
            numb=as.vector(tab)
            for (i in 1:length(loc)){
              lines(c(loc[i],loc[i]),c(b1,b1+step2*numb[i]),col='darkgrey')
            }
            points(loc,b1+step2*numb,col=2,pch=19)
            #text(loc,b1+step2*numb,loc,adj=0.5,cex=0.2)
          }
          
          
          if (length(gnohigh)>0){
            tab=table(data2[gnohigh,18])
            loc=as.numeric(names(tab))
            numb=as.vector(tab)
            for (i in 1:length(loc)){
              lines(c(loc[i],loc[i]),c(b2,b2-step2*numb[i]),col='darkgrey')
            }
            points(loc,b2-step2*numb,col='darkviolet',pch=19)
            #text(loc,b2-step2*numb,loc,adj=0.5,cex=0.2)
          }
          
          if (length(BLB)>0){
            tab=table(data2[BLB,18])
            loc=as.numeric(names(tab))
            numb=as.vector(tab)
            for (i in 1:length(loc)){
              lines(c(loc[i],loc[i]),c(b2,b2-step2*numb[i]),col='darkgrey')
            }
            points(loc,b2-step2*numb,col='darkgreen',pch=19)
            #text(loc,b2-step2*numb,loc,adj=0.5,cex=0.2)
          }
          
        }
        
        
        # # plot 4
        {
          data=databoth[which(databoth[,1]=="VUS" | databoth[,1]=="CON" ),]
          if(dim(data)[1]==0){
            text(max(maximum,1)/2,2.5+0.2,"No missense VUS and conflicting variant reported",adj=0.5,cex=1.5,col=2,font=3)
            text(max(maximum,1)/2,2.5-0.2,"No LoF VUS and conflicting variant reported",adj=0.5,cex=1.5,col=2,font=3)
          }
          if(dim(data)[1]>0){
            b1=2.5+demigene
            b2=2.5-demigene
            VUSm=which(data[,1]=="VUS" & data[,4]=="nonsynonymous SNV")
            CONm=which(data[,1]=="CON" & data[,4]=="nonsynonymous SNV")
            VUSl=which(data[,1]=="VUS" & (data[,4]=="frameshift deletion" | data[,4]=="frameshift insertion" | data[,4]=="stopgain"))
            CONl=which(data[,1]=="CON" & (data[,4]=="frameshift deletion" | data[,4]=="frameshift insertion" | data[,4]=="stopgain"))
            
            if((length(VUSm)+length(CONm))==0){
              text(max(maximum,1)/2,2.5+0.2,"No missense VUS and conflicting variant reported",adj=0.5,cex=1.5,col=2,font=3)
            }
            if((length(VUSl)+length(CONl))==0){
              text(max(maximum,1)/2,2.5-0.2,"No LoF VUS and conflicting variant reported",adj=0.5,cex=1.5,col=2,font=3)
            }
            
            lev=sort(unique(c(data[VUSm,18],data[CONm,18])))
            
            if(length(lev)>0){
              tab=table(factor(data[VUSm,18],levels=lev))
              loc=as.numeric(names(tab))
              numb=as.vector(tab)
              tab2=table(factor(data[CONm,18],levels=lev))
              numb2=as.vector(tab2)
              
              for (i in 1:length(loc)){
                if(numb[i]>0){
                  lines(c(loc[i],loc[i]),c(b1,b1+step4*numb[i]),col='darkgrey')
                  points(loc[i],b1+step4*numb[i],col=4,pch=19)
                }
                if(numb2[i]>0){
                  lines(c(loc[i],loc[i]),c(b1,b1+step4*numb2[i]),col='darkgrey')
                  points(loc[i],b1+step4*numb2[i],col=3,pch=19)
                }
              }
              #text(loc,b1+step4*numb,loc,adj=0.5,cex=0.2)
              #text(loc,b1+step4*numb2,loc,adj=0.5,cex=0.2)
            }
            
            lev=sort(unique(c(data[VUSl,18],data[CONl,18])))
            
            if(length(lev)>0){
              tab=table(factor(data[VUSl,18],levels=lev))
              loc=as.numeric(names(tab))
              numb=as.vector(tab)
              tab2=table(factor(data[CONl,18],levels=lev))
              numb2=as.vector(tab2)
              
              for (i in 1:length(loc)){
                if(numb[i]>0){
                  lines(c(loc[i],loc[i]),c(b2,b2-step4*numb[i]),col='darkgrey')
                  points(loc[i],b2-step4*numb[i],col='orange',pch=19)
                }
                if(numb2[i]>0){
                  lines(c(loc[i],loc[i]),c(b2,b2-step4*numb2[i]),col='darkgrey')
                  points(loc[i],b2-step4*numb2[i],col=2,pch=19)
                }
              }
              #text(loc,b2-step4*numb,loc,adj=0.5,cex=0.2)
              #text(loc,b2-step4*numb2,loc,adj=0.5,cex=0.2)
            }
            
          }
        }
        
        
        # plot 3
        {
          if(dim(cons)[1]>50){
            posi=1:length(consw)
            rem=which(posi>=zoom[1] & posi<=zoom[2])
            posi=(posi-zoom[1])/(zoom[2]-zoom[1])*maximum
            lines(posi[rem],consw[rem],col=3,lwd=3)
          }
          
          # PROTEIN DOMAINS
          
          # take id of protein
          {
            defaultW <- getOption("warn") 
            options(warn = -1) 
            id1=strsplit(as.character(isoform), split="\\.")[[1]][1]
            options(warn = defaultW)
            id2=refseq[which(refseq[,3]==id1),1]
            if (length(id2)>1){
              id3=""
              for (iso in 1:length(id2)){
                id3=c(id3,as.character(uniprot[which(as.character(uniprot[,1])==id2[iso]),2]))
              }
              if(length(which(id3==""))>0){
                id3=id3[-which(id3=="")]
                id3=id3[1]
              }
            } else {
              id3=as.character(uniprot[which(as.character(uniprot[,1])==id2),2])
            }
          }
          
          #id3=""
          
          # plot domains and region
          if(length(id3)>0){if(id3!="" & is.na(id3)=="FALSE"){
            
            # download data about the protein
            print(id3)
            rel_json <- drawProteins::get_features(id3)
            prot <- drawProteins::feature_to_dataframe(rel_json)
            prot=prot[sort(prot$begin, index.return=TRUE)$ix,]
            
            # computing size of annotations
            l=1:length(prot$type)
            for (g in 1:length(l)){
              l[g]=prot$end[g]-prot$begin[g]
            }
            
            # for domains and regions
            
            all=which(l>0 & (prot$type=="REGION" | prot$type=="DOMAIN" | prot$type=="ZN_FING" | prot$type=="REPEAT" | prot$type=="CHAIN" | prot$type=="TRANSIT" | prot$type=="DNA_BIND" | prot$type=="TRANSMEM" | prot$type=="TOPO_DOM") )
            prot[all,2]=gsub(" \\d+","",prot[all,2])
            prot[,2]=sapply(strsplit(as.character(prot[,2]), ";"), `[`, 1)
            if(length(unique(prot[all,2]))>=21){
              all=all[1:20] 
            }
            coldom=GetColors(n = length(unique(prot[all,2])), scheme = "discrete rainbow")
            
            
            
            # CHAINS
            dom1=which(l>0 &  prot$type=="CHAIN")
            if(length(dom1)>0){
              uni=unique(prot$description[dom1])
              ndom1=length(unique(prot$description[dom1]))
              cols=coldom[1:ndom1]
              for (d in 1:length(uni)){
                a=which(prot$description==uni[d])
                for (e in 1:length(a)){
                  prot$description[a[e]]="Mature chain"
                  
                  pos1=min(prot$begin[a[e]],max(maximum,1))
                  pos2=min(prot$end[a[e]],max(maximum,1))
                  if(pos1<zoom[1]){pos1=zoom[1]}
                  if(pos2>zoom[2]){pos2=zoom[2]}
                  
                  if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
                    pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
                    pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
                    rect(pos1,    3.5-demigene,   pos2,   3.5+demigene,   col=cols[d])
                  }
                }
              }
            } else {
              ndom1=0
            }
            
            # domains
            dom2=which(l>0 & (prot$type=="DOMAIN" | prot$type=="ZN_FING" | prot$type=="DNA_BIND"))
            prot[which(prot$type=="ZN_FING"),2]=paste("ZN_FING ",prot[which(prot$type=="ZN_FING"),2],sep="")
            if(length(dom2)>0){
              uni=unique(prot$description[dom2])
              ndom2=length(unique(prot$description[dom2]))
              cols=coldom[(ndom1+1):(ndom1+ndom2)]
              for (d in 1:length(uni)){
                a=which(prot$description==uni[d])
                for (e in 1:length(a)){
                  pos1=min(max(maximum,1),prot$begin[a[e]],na.rm=T)
                  pos2=min(prot$end[a[e]],max(maximum,1),na.rm=T)
                  if(pos1<zoom[1]){pos1=zoom[1]}
                  if(pos2>zoom[2]){pos2=zoom[2]}
                  
                  if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
                    pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
                    pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
                    rect(pos1,   3.5-demigene+0.015,  pos2 ,   3.5+demigene-0.015,col=cols[d])
                  }
                }
              }
              
            } else {
              ndom2=0
            }
            
            ### FOR REGIONs
            dom3=which(l>0 & (prot$type=="REGION"| prot$type=="REPEAT" | prot$type=="TRANSIT" | prot$type=="TRANSMEM" | prot$type=="TOPO_DOM"))
            if(length(dom3)>0){
              uni=unique(prot$description[dom3])
              ndom3=length(unique(prot$description[dom3]))
              cols=coldom[(ndom1+ndom2+1):(ndom1+ndom2+ndom3)]
              for (d in 1:length(uni)){
                a=which(prot$description==uni[d])
                for (e in 1:length(a)){
                  
                  pos1=min(prot$begin[a[e]],max(maximum,1))
                  pos2=min(prot$end[a[e]],max(maximum,1))
                  if(pos1<zoom[1]){pos1=zoom[1]}
                  if(pos2>zoom[2]){pos2=zoom[2]}
                  
                  if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
                    pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
                    pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
                    rect(pos1,   3.5-demigene+0.03,  pos2  ,   3.5+demigene-0.03,col=cols[d])
                  }
                }
              }
            }
            if(length(all)>0){
              text(0,3.35,"Mature chains (large), domains (medium) and regions (narrow) from Uniprot (>=10aa)",adj=0,cex=1.2)
              legend(0,3.28,unique(prot$description[c(dom1,dom2,dom3)]),fill=coldom,col=coldom,cex=0.5,ncol=3,bg='white')
            }
          }}
        }
        
        
        # clusters
        {
          clustPLP=0
          clustBLB=0
          
          clup=regionsp[which(regionsp[,2]==isoform),]
          if(length(which(regionsp[,2]==isoform))>1){
            clustPLP=1
            for (i in 1:dim(clup)[1]){
              pos1=as.numeric(clup[i,3])
              pos2=as.numeric(clup[i,4])
              if(pos1<zoom[1]){pos1=zoom[1]}
              if(pos2>zoom[2]){pos2=zoom[2]}
              if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
                pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
                pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
                rect(pos1,5.5+demigene-0.015,pos2,5.5+demigene+0.015,col='orange')
              }
              if(pos1==zoom[1] & pos2==zoom[2]){
                rect(0,5.5+demigene-0.015,maximum,5.5+demigene+0.015,col='orange')
              }
            }
          }
          if(length(which(regionsp[,2]==isoform))==1){
            clustPLP=1
            pos1=as.numeric(clup[3])
            pos2=as.numeric(clup[4])
            if(pos1<zoom[1]){pos1=zoom[1]}
            if(pos2>zoom[2]){pos2=zoom[2]}
            if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
              pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
              pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
              rect(pos1,5.5+demigene-0.015,pos2,5.5+demigene+0.015,col='orange')
            }
            if(pos1==zoom[1] & pos2==zoom[2]){
              rect(0,5.5+demigene-0.015,maximum,5.5+demigene+0.015,col='orange')
            }
          }
          
          club=regionsb[which(regionsb[,2]==isoform),]
          if(length(which(regionsb[,2]==isoform))>1){
            clustBLB=1
            for (i in 1:dim(club)[1]){
              pos1=as.numeric(club[i,3])
              pos2=as.numeric(club[i,4])
              if(pos1<zoom[1]){pos1=zoom[1]}
              if(pos2>zoom[2]){pos2=zoom[2]}
              if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
                pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
                pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
                rect(pos1,5.5-demigene-0.015,pos2,5.5-demigene+0.015,col='turquoise3')
              }
              if(pos1==zoom[1] & pos2==zoom[2]){
                rect(0,5.5-demigene-0.015,maximum,5.5-demigene+0.015,col='turquoise3')
              }
            }
          }
          if(length(which(regionsb[,2]==isoform))==1){
            clustBLB=1
            pos1=as.numeric(club[3])
            pos2=as.numeric(club[4])
            if(pos1<zoom[1]){pos1=zoom[1]}
            if(pos2>zoom[2]){pos2=zoom[2]}
            if((pos1>zoom[1] & pos1<zoom[2]) || (pos2>zoom[1] & pos2<zoom[2])){
              pos1=(pos1-zoom[1])/(zoom[2]-zoom[1])*maximum
              pos2=(pos2-zoom[1])/(zoom[2]-zoom[1])*maximum
              rect(pos1,5.5-demigene-0.015,pos2,5.5-demigene+0.015,col='turquoise3')
            }
            if(pos1==zoom[1] & pos2==zoom[2]){
              rect(0,5.5-demigene-0.015,maximum,5.5-demigene+0.015,col='turquoise3')
            }
          }
          
          # legend of clusters
          
          demiclust=0.015
          end=maximum/35
          beg=maximum/29
          
          if(clustBLB==1 & clustPLP==1){
            rect(0,5.05-demiclust,end,5.05+demiclust,col='turquoise3')
            rect(0,5.11-demiclust,end,5.11+demiclust,col='orange')
            text(beg,5.11,"Significant PLP clustering",adj=0,cex=0.9)
            text(beg,5.05,"Significant BLB clustering",adj=0,cex=0.9)
          }
          if(clustBLB==1 & clustPLP==0){
            rect(0,5.075-demiclust,end,5.075+demiclust,col='turquoise3')
            text(beg,5.075,"Significant BLB clustering",adj=0,cex=0.9)
          }
          if(clustBLB==0 & clustPLP==1){
            rect(0,5.075-demiclust,end,5.075+demiclust,col='orange')
            text(beg,5.075,"Significant PLP clustering",adj=0,cex=0.9)
          }
          if(clustBLB==0 & clustPLP==0){
            text(0,5.07,"No significant clustering detected",adj=0,cex=0.9)
          }
          
        }
        
      }
      
      time2=Sys.time()
      print(time2-time1)
      
      dev.off()
    }
  )
  
  output$downloadPlot2 <- downloadHandler(
    filename = function() { 
      
      isoform=input$isoform
      
      cchange=input$cpos
      
      pchange=input$ppos
      
      radio=input$radio
      
      if(pchange!="-" & radio=="Protein"){
        cchange=mut1[which(mut1[,1]==isoform & mut1[,4]==pchange)[1],2]
      }
      
      if(pchange=="-" & radio=="Protein"){
        cchange="-"
      }
      
      if(cchange!="-" & radio=="DNA"){
        pchange=mut1[which(mut1[,1]==isoform & mut1[,2]==cchange),4]
      }
      
      if(cchange=="-" & radio=="DNA"){
        pchange="-"
      }
      
      if(input$cpos!="-") {return(paste("MutScore_",input$gene,"_",input$isoform,"_",cchange,"_",pchange,'.pdf', sep=''))}
      if(input$cpos=="-") {return(paste("MutScore_",input$gene,"_",input$isoform,'.pdf', sep=''))}
    },
    content = function(file){
      pdf(file, height=8,width=14)
      
      print("Doing plot2 for PDF")
      time1=Sys.time()
      
      dataallgenes=variables$dataallgenes2
      bla=dataInput()
      isoform=bla[1]
      data=dataallgenes[which(dataallgenes[,2]==isoform),]
      
      if(dim(data)[1]>0){
        
        mut1=variables$mut2
        
        bla=dataInput()
        
        isoform=bla[1]
        zoom=c(as.numeric(bla[2]),as.numeric(bla[3]))
        constrack=bla[4]
        conswindow=as.numeric(bla[5])
        cchange=bla[6]
        pchange=bla[7]
        radio=bla[8]
        
        if(is.na(pchange)){pchange="-"}
        if(is.na(cchange)){cchange="-"}
        if(is.na(radio)){radio="DNA"}
        
        if(pchange!="-" & radio=="Protein"){
          cchange=mut1[which(mut1[,1]==isoform & mut1[,4]==pchange)[1],2]
        }
        
        if(pchange=="-" & radio=="Protein"){
          cchange="-"
        }
        
        if(cchange!="-" & radio=="DNA"){
          pchange=mut1[which(mut1[,1]==isoform & mut1[,2]==cchange),4]
        }
        
        if(cchange=="-" & radio=="DNA"){
          pchange="-"
        }
        
        load(file="PLP-BLB.RData")
        
        par(mfrow=c(3,1))
        hist(PLP,n=50,xlab="MutScore",cex.main=1.7,cex.lab=1.5,cex.axis=1.3,col='red',main="MutScore distribution for PLP missense variants from ClinVar",yaxt='n',ylab="")
        title(ylab="Frequency", line=(2), cex.lab=1.5)
        axis(2,at=c(0,3000,6000,9000,12000),labels=c(0,3000,6000,9000,12000),las=2,cex.axis=1.3,pos=0)
        
        if(cchange!="-" & pchange!="-"){
          muto=mut1[which(mut1[,2]==cchange & mut1[,1]==isoform),6]
          abline(v=muto,lwd=2,col='green')
          legend(0.5, 13000, legend=paste(dataallgenes[which(dataallgenes[,2]==isoform)[1],5],", ",isoform,":",cchange,", ",pchange,"    MutScore = ",muto,sep=""),col=c("green"), lty=1, cex=1.5, xjust=0.5,lwd=2,bg="white")
          
          n=50
          a=hist(PLP,n=n,plot=F)
          b=hist(BLB,n=n,plot=F)
          m=which.min(abs(a$mids-muto))[1]
          proba=round(a$density[m]/(a$density[m]+b$density[m])*100,digits=1)
          text(0.5,6000,paste("",proba,"% of variants with similar scores are PLPs",sep=""),cex=1.3)
          proba2=round(length(which(PLP>muto))/length(PLP)*100,digits=1)
          text(0.5,3500,paste("",proba2,"% of PLPs have a score higher than ",muto,sep=""),cex=1.3)
        }
        
        hist(BLB,n=50,xlab="MutScore",cex.main=1.7,cex.lab=1.5,cex.axis=1.3,col='blue',main="MutScore distribution for BLB missense variants from ClinVar",yaxt='n',ylab="")
        title(ylab="Frequency", line=(2), cex.lab=1.5)
        axis(2,at=c(0,1000,2000,3000),labels=c(0,1000,2000,3000),las=2,cex.axis=1.3,pos=0)
        if(cchange!="-" & pchange!="-"){
          abline(v=muto,lwd=2,col='green')
          legend(0.5, 3900, legend=paste(dataallgenes[which(dataallgenes[,2]==isoform)[1],5],", ",isoform,":",cchange,", ",pchange,"    MutScore = ",muto,sep=""),col=c("green"), lty=1, cex=1.5, xjust=0.5,lwd=2,bg="white")
          text(0.5,1800,paste("",100-proba,"% of variants with similar scores are BLBs",sep=""),cex=1.3)
          proba2=round(length(which(BLB<muto))/length(BLB)*100,digits=1)
          text(0.5,1100,paste("",proba2,"% of BLBs have a score lower than ",muto,sep=""),cex=1.3)
        }
        
      }
      
      time2=Sys.time()
      print(time2-time1)
      
      dev.off()
    }
  )
  
  
}

shinyApp(ui = ui, server = server)

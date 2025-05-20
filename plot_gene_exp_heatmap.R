library(DESeq2)
library(pheatmap)
library(tibble)
library(dplyr)

# Define a function to generate heatmap
plot_gene_expression_heatmap <- function(dds_deseq_all, gene_list, metadata, intgroup = "Treatment", scale_data = TRUE) {
  
  # Extract normalized counts
  norm_counts <- counts(dds_deseq_all, normalized = TRUE)
  
  # Subset only the selected genes
  gene_counts <- norm_counts[rownames(norm_counts) %in% gene_list, ]
  
  # Log-transform for better visualization
  gene_counts <- log2(gene_counts + 1) 
  
  # Convert to dataframe and add sample metadata
  gene_counts_df <- as.data.frame((gene_counts))  # Transpose for pheatmap
  
  # Remove metadata columns before plotting
#  gene_counts_df<-column_to_rownames(gene_counts_df,"Sample")
#  gene_counts_matrix <- as.matrix(gene_counts_df[, gene_list])
  
  # Generate heatmap
  pheatmap(gene_counts_df,
           annotation_colors = annotation_colors2,annotation_col = df2,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("#d4d000", "white", "#7603a3"))(50),
           show_rownames = TRUE,
           show_colnames = FALSE,
           main = "Moderate Gene Expression in anti-p40",
           scale="row")
  # group<-metadata$Treatment
  # col_anno <- HeatmapAnnotation(df = group,
  #                               col = list(Treatment = c("#02adba","#b52367","grey")))
  # col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  # Heatmap(gene_counts_df,
  #         row_names_gp = gpar(fontface = "italic")) 
  
}

#create spearet gene lists based on function. Include Saa1, Saa3
# Example usage
library(pheatmap)
retinol<-c("Retsat","Sult1a1","Rxrb","Abcg8","Abcg5","Rbp2","Rdh10","Bco2","Bco1","Rxrg","Dhrs3")
Lipid<-c("Dgkg", "Gal3st1", "Acss2", "Plppr2", "Pitpna", "Nr1i2", "Pla2g3", "Lrp2", "Aldh3b2", "Fads6", "Awat2", "Agmo", "Ttpa", "Acp6", "Hsd17b2", "Them4", "Rictor", "Ubl4a", "Cpt1a", "Phyh", "Acbd5", "Mecer", "Chka", "Pnpla8", "Cyp4b1", "Plaat3", "Ocrl", "Acsl3", "Hspg2", "Pm20d1", "Per2", "AclY", "Nceh1", "Agps", "Angptl8", "Dok7", "C1d", "Pccb", "Cyp1a1", "Pde3a", "Echdc3", "Chpt1", "Acot6", "Ppara", "Lpin2", "Plcd1", "B4galt4"
)
infl<-c("Gba1", "Ifi35", "Park7", "Metrnl", "Ets1", "Tnf", "Cdh5", "LgalS1", "Ccn4", "Il12b", "Sirpa", "Ccn3", "Vps35", "Map3k8", "Enpp3", "Ggt1", "Zbp1", "Acod1", "S100a9", "S100a8", "Ptges", "Sema7a", "Akna", "Dhx9", "Alox15", "Mgst2", "Lpl", "Fpr2", "Gata3", "Adamts12", "Mefv", "Pla2g7", "Fxr1", "Sbno2", "Sharpin", "NlRp6", "NlRp3", "Ripk1", "Apoe", "Snx6", "Gsdmd", "TgfB1", "NkG7", "Nmi", "Nfkb1", "Nfkbia", "CylD", "Lacc1", "Bcl6", "Fbxl2", "Myd88", "Il21", "Dagla", "Tnfaip3", "C1qtnf3", "Ufl1", "Casp12", "Mdk", "Hyal2", "Tmsb4x", "Casp4", "Casp1", "Ffar2", "Tnfaip8l2", "Jak2", "Trim65", "Gpsm3", "Ccr2", "Il10", "Anxa1", "Sphk1", "Mmp3", "Il16", "Tnfrsf1b", "Il17ra", "Tnfrsf1a", "Hck", "Psma6", "Ifng", "Il1b", "Tlr9", "Xcl1", "Tlr7", "Cd47", "Daglb", "Smpdl3b", "Dnase1l3", "Tlr4", "Tlr3", "BirC2", "Tlr2", "Vamp3", "Cebpa", "Grn", "Mvk", "Clec12a", "Gps2", "Nlrc3", "Aoah", "Nlrc4", "Nod2", "Ptgs2", "Usp18", "Rela", "Lilra5", "Mapk7", "Camk2n1", "Ccl5", "Fem1c", "Il33", "Il10ra", "Nr1h2", "Stat3", "Pla2g2a", "Nr1h4", "Nr1h3", "Mapk13", "NlRp10", "Bcl6b", "Otulin", "Ccm2", "Pik3ap1", "Ptpn2"
)
module1<-c("St6galnac3","St3gal1","St6galnac4","St3gal2","St6galnac5","Crtap","Adamts3","Bmp1","Col1a2","Col12a1","Col5a2","Loxl4","Loxl2","Ddr2","Mmp16","Ccdc80","Mmp2")
module2<-c("Dtna","Kcdn3","Tpm2","Tpm1","Lmod1","Tacr2","Slc8a1","Cacna1h","Myom2","Tpcn2","Mylk","Smtn","Des","Sgca","Asph","Fxyd1","Sspn","Myh11","Nmur1","Apbb1","Utrn","Cryab","Cda","Btg3","Btg2","Spin4","Tgfb1i1", "Jade1", "Mia3", "Spint2", "Rerg", "Sox17", "Csrp1", "Creb3l4", "Rps6ka2", "Bcl7c", "Dact3", "Speg", "Tns2", "Klf11", "Myocd", "Tesc", "Magi2", "Prox1", "Cggrf1", "Tent5b", "Atg13", "Terf2", "Cbfa2t3", "St7l", "Nupr1", "Cryab", "Plpp1", "Kank2", "Adam22", "Twist1", "Foxo4", "Dll1", "Thbs1", "Dnajb2", "Pdzd2", "Scin", "Ccdc85b", "Podxl", "Alox5", "Nacc2", "Smyd2", "Igfbp6", "Slit2", "Cgref1", "Mpp3", "Tgfb3", "Fuz", "Ptpn14", "Smarca2", "Mt3", "Tbx2", "Bmp5", "Bmp4", "Nr4a1", "Cdk6", "Proc", "Tmem98", "Tmem115", "Dlc")
module4<-c("Chrm2", "Sema5a", "Cntfr", "Nrsn1", "Myt1l", "Zmynd8", "Chd7", "Atn1", "Msi1", "Msi2", "Fgf1", "Vsig10", "Rapgefl1", "Ptprf", "Ophn1", "Dpysl5", "Dner", "Tnr", "Kat8", "Ephb1", "Mnx1", "Sh3gl2", "Cdon", "Rbfox1", "Sema6a", "Rbfox3", "Utp11", "Aplp2", "Pax6", "Nrg1", "Gfra1", "Nav2", "Gfra3", "Myt1", "Pou3f3", "Tgfbr1", "Thoc6", "Vcan", "Aldh5a1", "Nr5a2", "Adgrb2", "Adgrb1", "Kcnq2", "Tyro3", "Dcx", "Shank3", "Hdac4", "Nlgn1", "Chrna3", "Ctf1", "Nrxn1", "Wnt8b", "Nrxn2", "Ascl1", "Acvr1b", "Nr2c2", "Npas2", "Neurod1", "Rxrb", "Atxn3", "Dhx30", "Acvr1c", "Erbb3", "Btbd3", "Spock2", "St8sia4", "Spock1", "Igf2bp2", "Pcdh1", "Slc25a23", "Mark4", "Rxrg", "Gpm6b", "Atoh1", "Plxna3", "Mbd5", "Egr2", "Fn1", "Wwp1", "Cables1", "Grhl2", "Phox2b", "Ahi1", "Rnf103", "Fgf14", "Glrb", "Fgf18", "Pxmp2", "Pcdhb5", "Unc119", "Serpini1", "Cntn4", "Fgf13", "Fgf11", "Rb1", "Cntnap2", "Gsk3b", "Ptpru", "Tenm4", "Dock7", "Stmn2", "Ptprk", "Stmn1", "Micall1", "Rab6b", "Epha5", "Etv1", "Pbx1", "Runx1", "Gprin3", "Gprin2", "Ugdh", "Vapa", "Tdp2", "Wnk1", "Gprin1", "Rab35", "Tbc1d24", "Ulk1", "Insm1", "B4galt6", "Mapre2", "Cdk5r1", "Gria1", "Gria2", "Snap25", "Kcnc4", "Htr4", "Cplx1", "Grm7", "Ptchd1", "Kif5a", "Kcnmb4", "Slc12a6", "Slc12a7", "Ntsr1", "Chrnb2", "Unc13b", "Syt5", "Unc13c", "Syt1", "Oprk1", "Htr3a", "Grin2b", "Syn2", "Syn1", "Cd200r1", "Rit2", "Sypl1", "Amph", "Mycbpap")
module3<-c("Tbx2","Gpr107","Tpd52l1","Rem1","Mx1","Hnrnpd","Tubgcp3","Vps9d1","Col6a1","Ckb","Fndc5","Cnn1","Sgca","Ergic1","Stk38l","Rin2","Lrp3","Ltbp1","Nfix","Vsig2","Esam","Bbc3","Slc39a13","Rhoc","Ipo4","Slc25a42","Nherf2","Scin","Grin2d","Lamb1","Map2k7","Bcam","Prkar2b","Cyth2","Glg1","Dusp3","Fosb","Rnf215","Man1a","Syce2","Hlf","Cavin1","Map3k20","Pde1c","9330159f19rik","Large1","Rab23","Ulk2","Dgkq","Nes","Bcan","Crbn","Igf1r","Adcy9","Irag1","Kit","Reep5","Dnase1","Hspb7","Eps15l1","Upk1a","Crip2","Cbfa2t3","Tek","Neurl1a","Sfmbt1","Slc4a3","Usp19","Pola1","Saal1","Tmco6","P3h4","Ddah2","Atp1a2","Stx1a","Hapln4","Cav1","Hmgcll1","Phox2a","Fhl2","Fut1","Gprc5b","Trpm5","Sardh","Tmem115","Nxf1","Prox1","Optc","Rbm25","Kdelr3","Sugp1","Fuz","Evi5","Evi5l","Scarf2","Ubxn11","Amotl1","Zbtb47","Card14","Bcas1","Clip3","Hps5","Dll1","Nol3","Gata2","Cacfd1","Ttbk1","Abo","Sqstm1","Castor2","Chdh","Sertad4","Timp2","Atad5","Rhbdl3","Mmp9","Jph2","Ptgis","Cyth3","Cyb5r3","Ints2","Med13l","Mafk","Stk4","Tbx3","Copz2","Zfyve27","Myh11","P4ha2","Rab3d",
"Ppp1r12c","Ahr","Aoc3","Ffar3","Utrn","Smpd2","Wasf1","Traf3ip2","Arfgef3","Gopc","Lama2","Ube2d1","Rhobtb1","Tcp11l2","Timp3","Ddx50","Tacr2","Slc16a7","Slc1a4","Ptprb","Kcnmb1","Meis1","Col6a2","Nopchap1","Aldh1l2","Stk10","Xpo1","Asb3","Wdpcp","Sgcd","Pdlim4","Tns3","Btg2","Adcy1","Smtn","Cfap36","Sap30l","Galnt10","Myocd","Tspan13","Agr2","Klhl29","Trappc12","Klf11","Kif3c","Sntg2","Cep112","Rab37","Atp2a3","Tekt1","Mxra7","Rabep1","Kif1c","Rflnb","Spag9","Ankrd40","Dlg4","Per1","Myh10","Ntn1","Usp43","Krt36","Plcd3","Gtf2a1","Ptpn21","Mthfd1","Ppp2r5e","Ppm1a","Pcnx1","Nrde2","Golga5","Unc79","Tshz3","Tgfb3","Vash1","Kif26a","Gli3","Id4","Aspn","Fam217a","Eci3","Rpp40","Gadd45g","Auh","Ror2","Pdlim7","Rasa1","Agtpbp1","Dapk1","Nkd2","Pcsk1","Rhobtb3","Lpcat1","Xrcc4","Jmy","Plk2","Trim23","Ppwd1","Fgf10","Slc4a7","Plpp1","Gpx8","Nr1d2","Camk2g","Bmp4","Samd4","Rnase4","Galnt15","Chat","Pspc1","Gdf10","Lats2","Naa16","Scara5","Dpysl2","Ebf2","Bora","Fhip2b","Hr","Dmtn","Gfra2","Rbm26","Gpr180","Gdnf","Dab2","Ttc33","4931414p19rik","Pdzd2","Efs","Npr3","Rai14","Slc25a32","Zfpm2","Matn2","Sybu","Has2","Tab1","Dscc1","Fam118a","Rapgef3","Pmm1","Pde1b","Ppp1r1a","Plaat1","Tmem44","Slc52a2","Mapk12","Cblb","Retnlb","Boc","Zbtb20","Pmm2","Dgcr8","St3gal6","Tmem45a","Top3b","Lrch3","Mylk","Adcy5","Pdia5","Tmem41a","Tra2b","Btg3","Son","Paxbp1","Enah","Nr4a1","Krt18","Igfbp6","Fabp2","Gtf2ird1","Fhl1","Ache","Slc17a9","Nfatc4","Celsr3","Pcbp4","Ogfod2","Tiam2","Rps6ka2","Slc22a1","Igf2r","Smoc2","Ubr2","Kif6","Myom1","Clip4","Lbh","Ndufaf7","Cacna1h","Dnase1l2","Hagh","Tmem204","Luc7l","Dusp1","Spdef","Zfp523","Clps","Svil","Zeb1","Thada","Dtna")

FA<-c("Ehhadh", "Glul", "Car2", "Ppara", "Acsm3", "Reep6", "Retsat", "Hao2", "Prdx6", "Eno3", "Eno2", "Cpt1a", "Fabp2", "Hibch", "Ugdh", "Cyp1a1", "Bphl", "Acadm", "Decr1", "Hmgcs2", "Hadh", "Gcdh", "Acaa2", "Hpgd", "Hmgcs1", "Pdha1", "Acads", "Adipor2", "Suclg2", "Acss1", "Mcee", "Aldh3a1", "Ech1", "Urod", "Gstz1", "Acsl1", "Acaa1a", "Adh7", "Xist", "Cbr1", "Aadat", "Acot2", "Cpt2", "Ltc4s", "Inmt", "Etfdh", "Hsd17b10", "Crat", "D2hgdh", "Tdo2", "Fmo1", "Fasn", "Acadvl", "Hsp90aa1", "Ptprg", "Hadhb", "Echs1", "Gpd1", "Gapdhs", "Uros", "Hsph1", "Cd36", "Hsd17b4", "Sdha", "Aoc3", "Sucla2", "Suclg1", "Kmt5a", "Pcbd1", "Hsdl2", "Nthl1", "Dhcr24", "Maoa", "Cryz", "Bckdhb", "Dld", "Gabarapl1", "Hmgcl", "Acot8", "Fh1", "Erp29", "Cyp4a10", "Hsd17b7", "Apex1", "Acox1", "Eci1", "Nbn", "Cbr3", "Odc1", "Bmpr1b", "Sms-ps", "Idh3b", "Aco2", "Cidea", "Metap1", "Car4", "S100a10", "Cpox", "Aqp7", "Aldh3a2", "Dlst", "Auh", "Slc22a5", "G0s2", "Idh3g", "Gad2", "Mlycd", "Alad", "Pdhb", "Aldh9a1", "Adsl", "Sdhc", "Trp53inp2", "Serinc1", "Mgll", "Pts", "Ephx1", "Cd1d2", "Idh1", "Mdh1", "Ncaph2", "Mix23", "Rdh1", "Rdh11", "Hccs", "Ywhah", "Rap1gds1", "Nsdhl", "Eci2", "Mdh2", "Acsl4", "Idi1", "Aldh1a1", "Cel", "Acadl", "Aldoa", "Vnn1", "Fabp1", "Gpd2", "Car6", "Ldha", "Hsd17b11", "Sdhd", "Acat2", "Acsl5", "Psme1", "Ube2l6", "Elovl5", "Lgals1", "Blvra", "Mif", "Grhpr", "Ostc", "H2az1", "Me1"
)
module3_p2<-c("Dnajc18","Proc","Atp6v1g2","Zfp521","Lama3","Dpysl3","Rab27b","Spire1","Sncaip","Ldlrad4","Grpel2","Prrc1","Pdgfrb","Zp1","Syt7","Stambpl1","Tnfrsf25","Cabp4","Smarca2","Vldlr","Lgals12","Plaat5","Rbp4","Sorbs1","Tm9sf3","Pcgf6","Gk","Trub1","Pnliprp2","Hdgfl3","Cwf19l1","Maged2","Tro","Dnajc14","Baiap2","Prim1","Avil","Crisp3","Atp23","Inpp5a","Muc2","Dock9","Cbx2","Klhl4","Podxl","Alox5","Jade1","Clec3b","Rassf3","Prkar1b","Gria4","Sox17","Mettl21a","Adam23","Nrp2","Raph1","Trak2","Col3a1","Ercc5","Dst","Dnajb2","Speg","Des","Pask","Relch","Cd55","C4bp","Pigr","Csrp1","Slc45a3","Atp2b4","Cdc42bpa","Mptx1","Spta1","Tada1","Dpt","Smyd2","Ptpn14","Mark1","Suv39h2","Rgs5","Enkur","Pter","Dnajc1","Lypd6b","Abi2","Ak1","Nr4a2","Galnt5")
module3_p3<-c("Acvr1","Upp2","Zeb2","Gpsm1","Nacc2","Lrrc26","Hnmt","Slc43a1","Hoxd8","Pamr1","Meis2","Mapk8ip1","Atg13","Dnaaf9","Gfra4","Fermt1","Ttl","Bfsp1","Pag1","Snx16","Map1lc3a","Rbm39","Slc2a10","Pik3ca","Pex5l","Pld1","Trpc3","Postn","P2ry1","Tm4sf1","Wwtr1","Wnt2b","Tshb","Tspan2","Hmgcs2","Prpf38b","Slc6a17","Creb3l4","Slc50a1","Cyp2u1","Pla2g12a","Lrat","Gucy1b1","Snx7","Alpk1","Dkk2","Ptgfr","Khdc4","Myoz2","Ccn1","Asph","Faxc","Clca1","Odf2l","Rragd","Epha7","Cfap206","Col15a1","Tex10","Kif12","Tnc","Astn2","Kdm4c","Mpdz","Spink4","Ccdc107","Tpm2","Npr2","Mllt3","St3gal3","B4galt2","Bend5","Orc1","Edn2","Ccdc30","Trit1","Pomgnt1","Cda","Pink1","Gnat3","Ccdc28b","Prkag2","Slc4a2","Clcn6","Dffb","Mib2","Dvl1")
module3_p3<-c("Afap1","Pcdh7","Evc","Cgref1","Pi4k2b","Sel1l3","Cckar","Cenpc1","Gfi1","Tgfbr3","Zfp326","Sparcl1","Tesc","Pf4","Rilpl1")
module3_present <- module3[module3 %in% rownames(dds_deseq_all)]

gene_list_p40_ibd<-c("Igha","Ighg1","Il18","Defb45")
gene_list_1 <- c("Gzma","Gzmb","Il17a","Tnf","Il12b","Tbx21","Il1b","Tlr2","Irf5","Cxcl1","Ccl5","Mmp3","Lcn2","Cd44","Chit3l1","Slc11a1","Tgfb1","Nod2","Il22","Reg3b","Reg3g")  # p40 moderate and similar with IBD
gene_list_2<-c("Nos2","Nox1","Cybb","Duox2","Sod3","Gpx2","Gpx3","Prdx1","Prdx4","Prdx5") #p40 moderate ROS/RNS
gene_list_1_2<-c("Gzma","Gzmb","Il17a","Tnf","Il12b","Tbx21","Il1b","Tlr2","Irf5","Cxcl1","Ccl5","Mmp3","Lcn2","Cd44","Chit3l1","Slc11a1","Tgfb1","Nod2","Il22","Reg3b","Reg3g","Nos2","Nox1","Cybb","Duox2","Sod3","Gpx2","Gpx3","Prdx1","Prdx4","Prdx5")

gene_list_p40_c<-c("Ifng","Il12b","Cdkn1a","Mmp12","Cxcl2","Glb1","Mptx1","Mptx2","Reg4","Cyp2c55","Tlr9","Stat1","Stat3","Ltf") #similar to p40 and VC
gene_list_4<-c("Sync","Retnlb",) #Similar in 
gene_list_p40<-c("Osm","Il13ra2","Trem1","Il23a","Ccl20","Fpr2")
gene_list_non_resp_Li<-c("Socs3", "Cd55", "Kdm5d", "Igfbp5", "Slc15a1", "Xpnpep2", "Hla-dqa2", "Hmgcs2", "Ddx3y", "Itgb2", "Cdkn2b", "Hla-dqa1","Cxcl5")
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1401733/full
gene_list_non_resp_Granot <- c("F3","Cxcl1","Fam20a","Raph1","Ccl2","Cxcl3","Cxcl2","Pf4v1","Id1","Hp")

gene_list_epi<-c("Tff3","Spdef","Lgr5","Cldn1","Cldn7","Cldn8","Cldn4","Ocln","Muc2","Ccnd1") #no sig between p40 and control: Cldn4, Cdn1 up in p40, cdn8: up in control, 
plot_gene_expression_heatmap(dds_deseq_all, infl, metadata)

c("#02adba","#b52367","grey" )c("#02adba","#b52367","grey" )c("#02adba","#b52367","grey" )
annotation_colors2 <- list(
  Treatment = c(
    "HhWT+anti-p40" = "#02adba",
    "HhWT+isotype_control_mAb"="#b52367",
    "No_treatment_controls"="grey"))

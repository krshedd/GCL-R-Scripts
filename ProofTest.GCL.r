ProofTest.GCL=function(sillyvec,loci,groupnames,groupvec,samplesize,prefix,dir,prprtnl=TRUE,type="BAYES",suffix="",nreps=40000,nchains=5,priorvec,initmat,thin=c(1,1,1),switches="F T F T F T F"){
#########################################################################################################################################################################################################
#
#  This function removes individuals from the baseline to create a test mixture.
#
#  "type" is a vector containing any or all of the following:c("BAYES","BUGS","SPAM") 
#    If type="BUGS", the arguments "suffix", "nreps", "nchains", and "initmat" are ignored.
#
# Written by Jim Jasper 6/16/2011
# Changed by Kyle Shedd 1/9/2015: added "switches" and "thin" argument to allow changing CreateControlFile.GCL within this function and allow for Individual Assignment, changes made to line 236 to reflect
#########################################################################################################################################################################################################

  pt_prefix=prefix

  possibletypes=c("BAYES","BUGS","SPAM")

  empty=sapply(type,function(tp){
                                  if(is.na(match(tp,possibletypes))){
                                    stop(paste("What is type '",tp,"', hoser?  Should be either 'BAYES','BUGS',or 'SPAM'!!!",sep=""))
                                  } 
        })


  typeIND=sapply(possibletypes,function(tp){!is.na(match(tp,type))})

  mytypes=possibletypes[typeIND]
  
  if(!file.exists(dir)){stop(paste(dir," does not exist, hoser!!!"))}
  
  if(sum(is.na(match(loci,LocusControl$locusnames)))){stop(paste("'",loci[is.na(match(loci,LocusControl$locusnames))],"' from argument 'loci' not found in 'LocusControl' object!!!",sep=""))}

  nalleles=LocusControl$nalleles[loci]

  ploidy=LocusControl$ploidy[loci]

  mysillyvec=paste(pt_prefix,sillyvec,sep=".")
             
  empty=sapply(mysillyvec,function(silly){
                                          assign(paste(silly,".gcl",sep=""),get(paste(sillyvec[match(silly,mysillyvec)],".gcl",sep=""),pos=1),pos=1)
                                        }) 
 
  G=length(samplesize)

  names(samplesize)=groupnames

  if(G!=max(groupvec)){stop("samplesize is incompatible with groupvec, hoser!!!")}    
  
  Gvec=(1:G)[samplesize>0]
  Gnames=groupnames[Gvec]
  names(Gvec)=Gnames

  mysillys=sapply(Gnames,function(gname){
                    mysillyvec[groupvec==Gvec[gname]]
             },simplify=FALSE)
          
  my.gcl=sapply(Gnames,function(gname){
             sapply(mysillys[[gname]],function(silly){
               get(paste(silly,".gcl",sep=""),pos=1);
             },simplify=FALSE)
           },simplify=FALSE)

  sillydim=sapply(Gnames,function(gname){
               sapply(mysillys[[gname]],function(silly){
                 dim(my.gcl[[gname]][[silly]]$scores)
               },simplify=FALSE);
             },simplify=FALSE)


  n=sapply(Gnames,function(gname){
        sapply(mysillys[[gname]],function(silly){
                if(length(sillydim[[gname]][[silly]])>2){
                  ans=sillydim[[gname]][[silly]][1];
                };
                if(length(sillydim[[gname]][[silly]])==2){
                  ans=1;
                };
                if(length(sillydim[[gname]][[silly]])<2){
                  ans=NULL;
                };
                ans;
        },simplify=FALSE);
    },simplify=FALSE)

 
  IDs=sapply(Gnames,function(gname){
         sapply(mysillys[[gname]],function(silly){
                   mygcl=my.gcl[[gname]][[silly]]
                   counts=mygcl$counts;
                   attributes=mygcl$attributes;
                   if(length(sillydim[[gname]][[silly]])>2){          
                     ids=as.character(dimnames(counts)[[1]]);
                   };
                   if(length(sillydim[[gname]][[silly]])==2){          
                     ids=as.character(attributes$FK_FISH_ID);
                   };
                   if(length(sillydim[[gname]][[silly]])<2){          
                     ids=NULL;
                   };
                   ids;
                 },simplify=FALSE);
               },simplify=FALSE)

  isnaIDS=sapply(Gnames,function(gname){
                mylist=sapply(mysillys[[gname]],function(silly){
                         ans=sum(is.na(IDs[[gname]][[silly]]));
                         if(ans){
                           stop(paste("ID contain 'NA's for silly '",silly,"', hoser!!!",sep=""));
                         };
                         if(!ans){
                           NULL;
                         };
                        NULL;
                       },simplify=FALSE);
                NULL
            },simplify=FALSE)


  nsamps=sapply(Gnames,function(gname){
             nn=n[[gname]]  
             if(prprtnl){
               prob=nn
             };
             if(!prprtnl){
               prob=rep(1,length(nn))
             };
             mytab=table(sample(x=mysillys[[gname]],size=samplesize[gname],prob=prob,replace=TRUE));
             IND=match(names(mytab),mysillys[[gname]]);
             newtab=array(NA,c(length(mysillys[[gname]])),dimnames=list(mysillys[[gname]]));
             newtab[IND]=mytab;
             newtab[is.na(newtab)]=0;
             newtab
         },simplify=FALSE)
                  
  IDstosample=sapply(Gnames,function(gname){
                  sapply(mysillys[[gname]],function(silly){
                    if(nsamps[[gname]][silly]){
                      ans=sample(x=IDs[[gname]][[silly]],size=nsamps[[gname]][silly],replace=FALSE)
                    };
                    if(!nsamps[[gname]][silly]){
                      ans=NULL
                    };
                    ans},simplify=FALSE);
              },simplify=FALSE)


  mixsillys=sapply(Gnames,function(gname){
               array(paste("mix.",mysillys[[gname]],sep=""),c(length(mysillys[[gname]])),dimnames=list(mysillys[[gname]]))
            },simplify=FALSE)

  IDstoremove=sapply(Gnames,function(gname){
                sapply(mysillys[[gname]],function(silly){
                         if(nsamps[[gname]][silly]==0){
                           ID=IDs[[gname]][[silly]]
                         };
                         if(nsamps[[gname]][[silly]]){
                           ID=IDs[[gname]][[silly]][-match(IDstosample[[gname]][[silly]],IDs[[gname]][[silly]])]
                         };
                         ID
                       },simplify=FALSE);
              },simplify=FALSE)

  nIDsIND=sapply(Gnames,function(gname){
            sapply(mysillys[[gname]],function(silly){
              ans= length(IDstoremove[[gname]][[silly]]) + as.numeric(nsamps[[gname]][silly]) == as.numeric(n[[gname]][silly]);
              if(!ans){
                stop(paste("ID lengths don't agree for silly '",silly,"', hoser!!!",sep=""));
              }; 
              if(ans){
                NULL;
              };
              NULL 
            },simplify=FALSE)
            NULL     
          },simplify=FALSE)             

  empty=sapply(Gnames,function(gname){
          sapply(mysillys[[gname]],function(silly){
            mymixsilly=paste(mixsillys[[gname]][[silly]],".gcl",sep="");
            exsts=exists(mymixsilly);
            if(exsts){
              suppressWarnings(rm(list=mymixsilly))
            };
            assign(mymixsilly,get(paste(silly,".gcl",sep=""),pos=1),pos=1)
          })
        })

  empty=sapply(Gnames,function(gname){
          sapply(mysillys[[gname]],function(silly){
            mix.silly=mixsillys[[gname]][[silly]];
            mxids=list(IDstoremove[[gname]][[silly]]);
            names(mxids)=mix.silly;
            mydim=length(mxids);
            if(length(mydim)){              
              RemoveIDs.GCL(silly=mix.silly,IDs=mxids);
            };
          });
        }) ##  this is for mixture

  empty=sapply(Gnames,function(gname){
          sapply(mysillys[[gname]],function(silly){
            if(length(IDstosample[[gname]][[silly]])){
              bsids=list(IDstosample[[gname]][[silly]]);
              names(bsids)=silly;
              RemoveIDs.GCL(silly,bsids)
            }
           })
          }) ## this is for baseline  

    mixdir=sapply(mytypes[mytypes!="BUGS"],function(tp){paste(dir,"\\",tp,".mixture",sep="")})
    empty=sapply(mixdir,function(dr){if(!file.exists(dr)){dir.create(dr);print(paste("Created '",dr,"', hoser!!!",sep=""))}})    
    basedir=sapply(mytypes[mytypes!="BUGS"],function(tp){paste(dir,"\\",tp,".baseline",sep="")})
    empty=sapply(basedir,function(dr){if(!file.exists(dr)){dir.create(dr);print(paste("Created '",dr,"', hoser!!!",sep=""))}})
    outputdir=sapply(mytypes[mytypes!="BUGS"],function(tp){paste(dir,"\\",tp,".output",sep="")})
    empty=sapply(outputdir,function(dr){if(!file.exists(dr)){dir.create(dr);print(paste("Created '",dr,"', hoser!!!",sep=""))}})     
    controldir=sapply(mytypes[mytypes!="BUGS"],function(tp){paste(dir,"\\",tp,".control",sep="")})  
    empty=sapply(controldir,function(dr){if(!file.exists(dr)){dir.create(dr);print(paste("Created '",dr,"', hoser!!!",sep=""))}})
    mixIDs=vector("list",length(as.character(unlist(mixsillys))))
    names(mixIDs)=as.character(unlist(mixsillys))
    for(gname in Gnames){
      for(silly in mixsillys[[gname]]){
        myscores=get(paste(silly,".gcl",sep=""),pos=1)$scores
        ids=as.character(IDstosample[[gname]][[mysillys[[gname]][match(silly,mixsillys[[gname]])] ]])
        mylen=length(ids)
        if(mylen){
          mixIDs[[silly]]=rep(NA,mylen)
          mixIDs[[silly]][1:mylen]=ids
        }        
      } 
    } 

    mixfortran=sapply(mytypes[mytypes!="BUGS"],function(tp){CreateMixture.GCL(sillys=as.character(unlist(mixsillys)),loci,IDs=mixIDs,mixname=pt_prefix,dir=mixdir[tp],type=tp)},simplify=FALSE)

    basefortran=sapply(mytypes[mytypes!="BUGS"],function(tp){CreateBaseline.GCL(sillyvec=mysillyvec,loci,dir=basedir[tp],basename=pt_prefix,type=tp,groupvec=groupvec)},simplify=FALSE)

    if(!is.na(match("BAYES",mytypes))){
      empty=CreateControlFile.GCL(sillyvec=mysillyvec,loci=loci,mixname=pt_prefix,basename=pt_prefix,suffix="",nreps=nreps,nchains=nchains,groupvec=groupvec,priorvec=priorvec,initmat=initmat,dir=controldir["BAYES"],thin=thin,mixfortran=mixfortran[["BAYES"]],basefortran=basefortran[["BAYES"]],switches=switches)
    }
    if(!is.na(match("SPAM",mytypes))){
      empty=CreateSPAMcontrol.GCL(sillyvec=mysillyvec,groupnames=groupnames,mixvec=pt_prefix,groupvec=groupvec,baselinefile=pt_prefix,loci=loci,dir=controldir["SPAM"],outdir=outputdir["SPAM"],mixdir=mixdir["SPAM"],basedir=basedir["SPAM"])
    }

  if(!is.na(match("BUGS",mytypes))){
    mixname="temp"  
    if(!exists("writeDatafileR")){source("V:\\DATA\\R_GEN\\nonGCL Source Scripts\\writeDatafileR.r")}
    BUGSdir=paste(dir,"\\BUGS",sep="")
    ToWhere=paste(BUGSdir,"\\",pt_prefix,"BUGSdata.txt",sep="")
    if(!file.exists(BUGSdir)){dir.create(BUGSdir);print(paste("Created '",BUGSdir,"', hoser!!!",sep=""))}    
    PoolCollections.GCL(collections=as.character(unlist(mixsillys)),loci=loci,IDs=NULL,newname=mixname)
    x=get(paste(mixname,".gcl",sep=""),pos=1)$counts[,loci,]
    y=FreqPop.GCL(mysillyvec,loci)
    n=t(sapply(mysillyvec,function(silly){sapply(loci,function(locus){sum(y[silly,locus,1:nalleles[locus]])})}))
    ORD=order(groupvec)
    Cg=as.vector(table(groupvec[ORD]))
    cmCg=c(0,cumsum(Cg))
    data=list(C=length(sillyvec),cmCg=cmCg,G=G,D=length(loci),Jd=as.numeric(nalleles),ploidy=as.numeric(ploidy),M=nrow(x),x=x,y=y[ORD,,],n=n[ORD,],alpha=as.numeric(priorvec[ORD]))
    writeDatafileR(DATA=data,towhere=ToWhere)
    suppressWarnings(rm(list=c(paste(mixname,".gcl",sep=""))))
  }


  suppressWarnings(rm(list=c(paste(mysillyvec,".gcl",sep=""),paste(unlist(mixsillys),".gcl",sep="")),pos=1))

  return(IDstosample)
}

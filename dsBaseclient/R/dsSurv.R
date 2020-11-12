
#'
#' @title Tests for correlation between paired samples
#' @description This function is similar to R function \code{cor.test}.
#' @details The function runs a two-sided Pearson test with a 0.95 confidence level.
#' @param x a character string providing  the name of a numerical vector.
#' @param y a character string providing  the name of a numerical vector.
#' @return the results of the survival object.
#' @author Sofack, Ghislain. (based on corTestDS by Demetris Avraam, for DataSHIELD Development Team)
#' @export
#'






ds.Surv <- function(x=NULL, y= NULL, newobj= NULL, datasources=NULL){

  # if no opal login details are provided look for 'opal' objects in the environment


  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  if(is.null(x)){
    stop("x=NULL. Please provide the names of the 1st numeric vector!", call.=FALSE)
  }
  if(is.null(y)){
    stop("y=NULL. Please provide the names of the 2nd numeric vector!", call.=FALSE)
  }




  # the input variable might be given as column table (i.e. D$object)
  # or just as a vector not attached to a table (i.e. object)
  # we have to make sure the function deals with each case

  objects <- c(x,y)
  xnames <- dsBaseClient:::extract(objects)
  varnames <- xnames$elements
  obj2lookfor <- xnames$holders

  # check if the input object(s) is(are) defined in all the studies
  for(i in 1:length(varnames)){
    if(is.na(obj2lookfor[i])){
      defined <- dsBaseClient:::isDefined(datasources, varnames[i])
    }else{
      defined <- dsBaseClient:::isDefined(datasources, obj2lookfor[i])
    }
  }



  # call the internal function that checks the input object(s) is(are) of the same class in all studies.
  for(i in 1:length(objects)){
    typ <- dsBaseClient:::checkClass(datasources, objects[i])
  }



  # create a name by default if user did not provide a name for the new variable
  if(is.null(newobj)){
    newobj <- "Survobj"
  }




  # call the server side function
  calltext <- call("SurvDS", x, y)

  DSI::datashield.assign(datasources, newobj, calltext)



  #############################################################################################################
  #DataSHIELD CLIENTSIDE MODULE: CHECK KEY DATA OBJECTS SUCCESSFULLY CREATED                                  #
  #
  #SET APPROPRIATE PARAMETERS FOR THIS PARTICULAR FUNCTION                                                 	#
  test.obj.name<-newobj
  #
  #
  # CALL SEVERSIDE FUNCTION                                                                                	#
  calltext <- call("testObjExistsDS", test.obj.name)
  #
  object.info<-DSI::datashield.aggregate(datasources, calltext)
  #
  # CHECK IN EACH SOURCE WHETHER OBJECT NAME EXISTS
  # AND WHETHER OBJECT PHYSICALLY EXISTS WITH A NON-NULL CLASS
  num.datasources<-length(object.info)
  #
  #
  obj.name.exists.in.all.sources<-TRUE
  obj.non.null.in.all.sources<-TRUE
  #
  for(j in 1:num.datasources){
    if(!object.info[[j]]$test.obj.exists){
      obj.name.exists.in.all.sources<-FALSE
    }
    if(is.null(object.info[[j]]$test.obj.class) || object.info[[j]]$test.obj.class=="ABSENT"){														 	#
      obj.non.null.in.all.sources<-FALSE
    }
  }

  if(obj.name.exists.in.all.sources && obj.non.null.in.all.sources){

    return.message<-
      paste0("A data object <", test.obj.name, "> has been created in all specified data sources")		 	#


  }else{

    return.message.1<-
      paste0("Error: A valid data object <", test.obj.name, "> does NOT exist in ALL specified data sources")	#

    return.message.2<-
      paste0("It is either ABSENT and/or has no valid content/class,see return.info above")				 	#

    return.message.3<-
      paste0("Please use ds.ls() to identify where missing")


    return.message<-list(return.message.1,return.message.2,return.message.3)

  }

  calltext <- call("messageDS", test.obj.name)
  studyside.message<-DSI::datashield.aggregate(datasources, calltext)

  no.errors<-TRUE
  for(nd in 1:num.datasources){
    if(studyside.message[[nd]]!="ALL OK: there are no studysideMessage(s) on this datasource"){			#
      no.errors<-FALSE
    }
  }


  if(no.errors){
    validity.check<-paste0("<",test.obj.name, "> appears valid in all sources")							    #
    return(list(is.object.created=return.message,validity.check=validity.check))						    #
  }

  if(!no.errors){
    validity.check<-paste0("<",test.obj.name,"> invalid in at least one source. See studyside.messages:")   #
    return(list(is.object.created=return.message,validity.check=validity.check,					    		#
                studyside.messages=studyside.message))			                                            #
  }


  #END OF CHECK OBJECT CREATED CORECTLY MODULE
  #############################################################################################################

}


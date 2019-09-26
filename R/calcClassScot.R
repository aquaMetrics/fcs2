#' Calculate Classification Results for Scotland
#'
#' @encoding UTF-8
#' @title Calculate Classification Results for Scotland
#'
#' @description
#' R Script for calculating new EQRs (for latest FCS2 model) for Scotland
#'
#' @param data dataframe with columns: \describe{ \item{DataHeldBy}{Data owner
#'   descriptor	Field data} \item{SiteCode}{Unique site code e.g. SEPA location
#'   code	Field data} \item{SiteName}{Site name	Field data}
#'   \item{SurveyDate}{Date of survey	Field data} \item{WBId}{SEPA Water body
#'   ID	SEPA GIS system/ SEWeb} \item{WBName}{SEPA Waterbody Name	SEPA GIS
#'   system/ SEWeb} \item{NumberOfRuns}{Number of separate runs in
#'   electrofishing survey	Field data} \item{SurveyArea}{Wetted area of survey
#'   site in square metres	Field data 	5-1628} \item{WetWidth}{Wetted width of
#'   survey site in metres	Field data 	0.3-41} \item{Slope}{Slope of survey
#'   site, in metres per kilometre	GIS or paper maps	0.01-153}
#'   \item{BarrierType}{Type of impassable barrier downstream, if present. 3
#'   possible values: Artificial, Manmade or Blank	Field data  and/or analysis
#'   of barrier data on SEWeb} \item{ImpassableBarriers}{Is the site above an
#'   impassable barrier? 2 possible values; 1= present, 0= absent	Field data
#'   and/or analysis of barrier data on SEWeb}
#'   \item{CatchmentAreaUpstream}{Total catchment area upstream in square
#'   km	GIS	0.05-2550} \item{CatchmentDrainageDirection}{Catchment drainage
#'   direction. Where does catchment enter the sea? Seven possible categories:
#'   SE= Tweed estuary to Fife Ness E=Fife Ness to Fraserburgh NE= Fraserburgh
#'   to Duncansby Head N= Duncansby Head to Cape Wrath NW= Cape Wrath to Mallaig
#'   W= Mallaig to Mull of Kintyre SW= Mull of Kintyre to Solway Firth}
#'   \item{GeologyClass}{Predominant geology type in catchment. 3 possible
#'   categories: S=Siliceous, C=Calcareous, O=Organic} \item{Altitude}{Site
#'   altitude (metres above sea level)	Field data  or GIS	1-499}
#'   \item{DistanceFromSource}{Distance from the site to the furthest upstream
#'   source (metres)	GIS	115-125066} \item{DistanceToSea}{Distance downstream
#'   from the site to the tidal limit (metres)	GIS	2-160267}
#'   \item{AnnualMeanFlow}{	Annual mean flow in cubic metres per
#'   second	Hydrological modelling (Low Flows Enterprise of similar)
#'   0.003-59.9} \item{AlkalinityValue}{	Alkalinity as CaCO3 (mg/L) 	Chemistry
#'   monitoring data	0.54-263.11} \item{TotalPValue}{	Total Phosphorus as P
#'   (mg/L) 	Chemistry monitoring data	0.006-0.505} \item{DOCValue}{	Dissolved
#'   organic carbon<0.45\eqn{\mu}m (mg/L)	Chemistry monitoring data	2.023-24.2}
#'   \item{SuspendedSolidsValue}{	Suspended solids (mg/L)	Chemistry monitoring
#'   data	0.63-32.36} \item{HydrometricAreaNo}{	The SEPA-defined hydrometric
#'   area which the site lies in	SEPA GIS layer}
#'   \item{LandUse.AgriculturalAreas}{% of catchment with Agricultural landuse.
#'   Defined as the following Corrine codes: Pastures (231), Non-irrigated
#'   arable land (211), Complex cultivation patterns (242), Land principally
#'   occupied by agriculture with significant areas or natural vegetation (243)
#'   GIS analysis of Corrine data	0-99.9} \item{LandUse.ConiferousForests}{% of
#'   catchment with coniferous forest landuse (Corrine Code 312)	GIS analysis of
#'   Corrine data	0-92.9} \item{Field}{Description	Source}
#'   \item{LandUse.Wetlands}{% of catchment with wetland landuse. Defined as the
#'   following Corrine codes: Inland marshes (411), Peat bogs (412), Intertidal
#'   flats (423)	GIS analysis of Corrine data	0-95.7}
#'   \item{Substrate.Small}{Total % of survey bed substrate composed of "High
#'   Organic", "silt", "sand" and "gravel", from SFCC survey protocol	Field data
#'   0-100} \item{Substrate.Large}{Total % of survey bed substrate composed of
#'   "Pebble", "cobble" and "boulder" from SFCC survey protocol	Field data
#'   0-100} \item{Substrate.Bedrock}{Total % of survey bed composed of "Bedrock"
#'   from SFCC survey protocol	Field data 	0-90}
#'   \item{Salmon_fry.Run1Total}{Total salmon fry (0+) caught in Run 1	Field
#'   data 	0-1014} \item{Salmon_fry.Run2Total}{Total salmon fry (0+) caught in
#'   Run 2. Leave blank if no 2nd run fished	Field data 	0-352}
#'   \item{Salmon_fry.Run3Total}{Total salmon fry (0+) caught in Run 3. Leave
#'   blank if no 3rd run fished	Field data 	0-251}
#'   \item{Salmon_fry.Run4Total}{Total salmon fry (0+) caught in Run 4. Leave
#'   blank if no 4th run fished	Field data 	0-12}
#'   \item{Salmon_parr.Run1Total}{Total salmon parr (1++) caught in Run 1	Field
#'   data 	0-133} \item{Salmon_parr.Run2Total}{Total salmon parr (1++) caught in
#'   Run 2. Leave blank if no 2nd run fished	Field data 	0-68}
#'   \item{Salmon_parr.Run3Total}{Total salmon parr (1++) caught in Run 3. Leave
#'   blank if no 3rd run fished	Field data 	0-34}
#'   \item{Salmon_parr.Run4Total}{Total salmon parr (1++) caught in Run 4. Leave
#'   blank if no 4th run fished	Field data 	0-15}
#'   \item{Trout_fry.Run1Total}{Total salmon fry (0+) caught in Run 1	Field data
#'   0-508} \item{Trout_fry.Run2Total}{Total salmon fry (0+) caught in Run 2.
#'   Leave blank if no 2nd run fished	Field data 	0-225}
#'   \item{Trout_fry.Run3Total}{Total salmon fry (0+) caught in Run 3. Leave
#'   blank if no 3rd run fished	Field data 	0-96}
#'   \item{Trout_fry.Run4Total}{Total salmon fry (0+) caught in Run 4. Leave
#'   blank if no 4th run fished	Field data 	0-59}
#'   \item{Trout_parr.Run1Total}{Total salmon parr (1++) caught in Run 1	Field
#'   data 	0-120} \item{Trout_parr.Run2Total}{Total salmon parr (1++) caught in
#'   Run 2. Leave blank if no 2nd run fished	Field data 	0-70}
#'   \item{Trout_parr.Run3Total}{Total salmon parr (1++) caught in Run 3. Leave
#'   blank if no 3rd run fished	Field data 	0-43}
#'   \item{Trout_parr.Run4Total}{Total salmon parr (1++) caught in Run 4. Leave
#'   blank if no 4th run fished	Field data 	0-25}
#'   }
#' @return dataframe
#' @export
#'
#' @examples
#'   \dontrun{
#'   results <- calcClassScot(data = fcs2::demo_data)#'
#' }
calcClassScot <- function(data) {
  colnames(data) <- sub("#", ".", colnames(data))

  ## Set pressure variables to reference values (values set using database version 13 - see section 6.4 of Phase 3 report)
  ##
  cat("Setting pressure variables to reference values\n")

  # Ammonium:
  subset <- !is.na(data$AlkalinityValue) & (data$AlkalinityValue < 50 |
                  (!is.na(data$Altitude) &  data$Altitude > 80 &
               data$AlkalinityValue < 200))

  data$AmmoniumValue[subset] <- 0.02362857

  subset <- !is.na(data$AlkalinityValue) & (data$AlkalinityValue > 200 |
                  (!is.na(data$Altitude) & data$Altitude <= 80 &
               data$AlkalinityValue > 50))

  data$AmmoniumValue[subset] <- 0.06622857

  subset <- is.na(data$AlkalinityValue) | (is.na(data$Altitude) &
           !is.na(data$AlkalinityValue) & data$AlkalinityValue > 50 &
              data$AlkalinityValue < 200)
  # kill Ammonium values where alk or alt unknown (when needed)
  data$AmmoniumValue[subset] <- NA

  # BOD:
  subset <- !is.na(data$AlkalinityValue) & (data$AlkalinityValue < 50 |
                  (!is.na(data$Altitude) & data$Altitude > 80 &
               data$AlkalinityValue < 200))

   data$BODValue[subset] <- 1.170833

  subset <- !is.na(data$AlkalinityValue) & (data$AlkalinityValue > 200 |
                  (!is.na(data$Altitude) & data$Altitude <= 80 &
               data$AlkalinityValue > 50))

  data$BODValue[subset] <- 1.645

  subset <- is.na(data$AlkalinityValue) | (is.na(data$Altitude) &
           !is.na(data$AlkalinityValue) & data$AlkalinityValue > 50 &
              data$AlkalinityValue < 200)
  # kill BOD values where alk or alt unknown (when needed)
  data$BODValue[subset] <- NA


  # ImpassableBarriers:  (only natural barriers at reference)
  subset <- !is.na(data$ImpassableBarriers) & !is.na(data$BarrierType) &
               data$ImpassableBarriers == 1 &
          (data$BarrierType == "Artificial" |
               data$BarrierType == "Manmade")

  data$ImpassableBarriers[subset] <- 0
  data$BarrierType[subset] <- NA

  subset <- !is.na(data$ImpassableBarriers) & !is.na(data$BarrierType) &
               data$ImpassableBarriers == 1 &
    (data$BarrierType == "Natural and manmade")

  data$BarrierType[subset] <- "Natural"

  subset <- is.na(data$ImpassableBarriers) | (is.na(data$BarrierType) &
           !is.na(data$ImpassableBarriers) &  data$ImpassableBarriers == 1)

  data$ImpassableBarriers[subset] <- NA
  data$BarrierType[subset] <- NA

  # LandUse.AgriculturalAreas:
  data$LandUse.AgriculturalAreas <- 0.1722859

  # LandUse.ConiferousForests:
  data$LandUse.ConiferousForests <- 0.58

  # SalmonStocking
  data$SalmonStocking <- 0 # no salmon stocking at reference

  # SRPMRP
  subset <- !is.na(data$AlkalinityValue) & !is.na(data$Altitude) &
               data$AlkalinityValue < 50 & data$Altitude <= 80

  data$SRPMRPValue[subset] <- 0.00971875

  subset <- !is.na(data$AlkalinityValue) & !is.na(data$Altitude) &
    data$AlkalinityValue < 50 & data$Altitude > 80

  data$SRPMRPValue[subset] <- 0.01025

  subset <- !is.na(data$AlkalinityValue) & data$AlkalinityValue > 50
  data$SRPMRPValue[subset] <- 0.03106667

  subset <- is.na(data$AlkalinityValue) | (!is.na(data$AlkalinityValue) &
              data$AlkalinityValue < 50 & is.na(data$Altitude))
  # kill SRPMRP values where alk or alt unknown (when needed)
  data$SRPMRPValue[subset] <- NA

  # pH
  data$pHValue <- 7.41

  ## Additional changes to data
  ##

  # set Slope and Distances < 1 to 1 so can log them
  data$Slope[data$Slope < 1] <- 1
  data$DistanceToSea[data$DistanceToSea < 1] <- 1
  data$DistanceFromSource[data$DistanceFromSource < 1] <- 1
  data$CatchmentAreaUpstream[data$CatchmentAreaUpstream < 1] <- 1

  # create additional column from WBId and SiteCode (used as a proxy for WBID
  # to combine surveys since WBId may be missing)
  joinByVar <- as.character(data$WBId)
  joinByVar[is.na(joinByVar)] <- as.character(data$SiteCode[is.na(joinByVar)])
  data$JoinByVar <- factor(joinByVar)


  ## Save data at reference (as R workspace)
  #cat("Saving data at reference\n")
  #save(file = "ScotDataAtRef.RData", data)


  ## Load all model fits
  ## WARNING: Loading all at once takes up a lot of memory!
  cat("Loading model fits\n")
  flush.console()

  SalmonFryFit <- fcs2:::SalmonFryFit
  SalmonParrFit <- fcs2:::SalmonParrFit
  TroutFryFit <- fcs2:::TroutFryFit
  TroutParrFit <- fcs2:::TroutParrFit
  ## Calculate all EQRs (single and joint EQRs for each survey and joined
  ## by waterbody)
  ## NOTE: This may take some time! - reduce n.samples and n.sims to speed
  ## up but lose accuracy
  ##
  cat("Calculating EQRs\n")
  flush.console()
  eqrs <- fcs2JointAndSingleEQR(SalmonFryFit,
                                SalmonParrFit,
                                TroutFryFit,
                                TroutParrFit,
                                newData = data,
                                joinByVar = "JoinByVar",
                                both = TRUE,
                                na.action = na.pass,
                                n.samples = 500,
                                n.sims = 1000,
                                showProgress = TRUE
  )

  ## Save EQRs with boundaries
  ##
  #cat("Saving EQRs\n")

  # extract EQR objects from list 'eqrs'
  EQRBySurvey <- eqrs$EQRBySurvey # EQRs by survey
  EQRByWBId <- eqrs$EQRByWBId # EQRs by waterbody

  ## NOTE: EQRs share order of database object 'data' so compare the two
  ## to identify surveys/sites

  # set EQR boundary
  boundaries <- c(0.009, 0.11, 0.6, 0.845) # calibrated boundaries for Scotland

  # save
  # save(file = "EQRs.RData", EQRBySurvey, EQRByWBId, boundaries)

  ## Summarise EQRs
  ##
  cat("Calculating EQR summary\n")
  flush.console()
  eqrSum <- fcs2EQRSummaryMatrix(SalmonFryFit,
                                 SalmonParrFit,
                                 TroutFryFit,
                                 TroutParrFit,
                                 newData = data,
                                 joinByVar = "JoinByVar",
                                 na.action = na.pass,
                                 predictions = "all",
                                 eqrs = eqrs,
                                 boundaries = boundaries
  )

  # correct column names
  # better describe class probs
  colnames(eqrSum)[2:6] <- paste("All species WB EQR",
                                 colnames(eqrSum)[2:6], "%")
  colnames(eqrSum) <- sub("JoinByVar",
                          "WB",
                          colnames(eqrSum)) # replace 'JoinByVar' with 'WB'

  # add WFD classification by survey
  classProbs <- t(fcs2Classify(EQRBySurvey[, rownames(eqrSum), 1],
                               boundaries = boundaries))
  colnames(classProbs) <- paste("All species survey EQR",
                                colnames(classProbs), "%")
  eqrSum <- cbind(
    eqrSum[, 1:7],
    classProbs,
    eqrSum[, -(1:7)]
  )

  # add data columns
  eqrSum <- cbind(data[rownames(eqrSum), c("DataHeldBy",
                                           "SurveyDate",
                                           "SiteCode",
                                           "SiteName",
                                           "WBId",
                                           "WBName")],
                  eqrSum[, -1])

  # Add row names and sort
  eqrSum <- eqrSum[match(1:nrow(eqrSum), rownames(eqrSum)), ]

  ## Save EQR summary
  ##
  cat("return EQR summary\n")

  # write to CSV
  ## NOTE: file should now be tidied up manually in Excel
  #write.csv(file = "EQRSummary.csv", eqrSum, row.names = TRUE, na = "")

  # save to R workspace
  #save(file = "EQRSummary.RData", eqrSum)
  #cat("Done\n")
  return(eqrSum)
}

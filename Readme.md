 # qPCR data
Copied columns A-L from "Phase III cold_8_12.xls" tab "Combined"
-deleted 16 rows at end of sheet that contained no S/H data
Copied columns A-H and O-Q from "Phase III hot_8_12.xls" tab "Combined"
Combined data into single sheet, added a column named "Temp" with value "cool" for corals from cooling treatment, and "heat" for corals from heating treatment
-filled in a few missing values for "mother" based on sample name
Wrote to file names "phaseIII_qPCR.csv"
-replaced 0.0001 and 0.00001 values for t6 S:H ratios with zeros, double checking that they should be zero based on qPCR Ct values. Rachel had previously replaced these zeros with small numbers to allow log transformation.
Further manipulations done in R, see script "data_carpentry.R"

# PAM data
Copied columns A-J and L-N from "Phase 3 composite.xls"

manually removed several rows with no data from "phaseIII_PAM.csv"


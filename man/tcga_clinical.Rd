\name{tcga_clinical}
\alias{tcga_clinical}
\docType{data}
\title{
  TCGA Clinical data
}
\description{
  TCGA Clinical data used on TCGA Explorer downloaded from GDC v32.0.
}
\usage{data("tcga_clinical")}
\format{
  A data frame with 10222 observations on the following 98 variables.
  \describe{
    \item{\code{barcode}}{a character vector}
    \item{\code{ajcc_pathologic_stage}}{a character vector}
    \item{\code{tissue_or_organ_of_origin}}{a character vector}
    \item{\code{status}}{a character vector}
    \item{\code{time}}{a numeric vector}
    \item{\code{days_to_last_follow_up}}{a numeric vector}
    \item{\code{age_at_diagnosis}}{a numeric vector}
    \item{\code{primary_diagnosis}}{a character vector}
    \item{\code{year_of_diagnosis}}{a numeric vector}
    \item{\code{ajcc_pathologic_t}}{a character vector}
    \item{\code{morphology}}{a character vector}
    \item{\code{ajcc_clinical_m}}{a character vector}
    \item{\code{ajcc_pathologic_n}}{a character vector}
    \item{\code{icd_10_code}}{a character vector}
    \item{\code{site_of_resection_or_biopsy}}{a character vector}
    \item{\code{tumor_grade}}{a character vector}
    \item{\code{cigarettes_per_day}}{a numeric vector}
    \item{\code{alcohol_history}}{a character vector}
    \item{\code{years_smoked}}{a numeric vector}
    \item{\code{alcohol_intensity}}{a logical vector}
    \item{\code{race}}{a character vector}
    \item{\code{gender}}{a character vector}
    \item{\code{ethnicity}}{a character vector}
    \item{\code{vital_status}}{a character vector}
    \item{\code{days_to_birth}}{a numeric vector}
    \item{\code{year_of_birth}}{a numeric vector}
    \item{\code{year_of_death}}{a numeric vector}
    \item{\code{days_to_death}}{a numeric vector}
    \item{\code{ajcc_pathologic_m}}{a character vector}
    \item{\code{irs_stage}}{a logical vector}
    \item{\code{iss_stage}}{a logical vector}
    \item{\code{ann_arbor_clinical_stage}}{a character vector}
    \item{\code{enneking_msts_stage}}{a logical vector}
    \item{\code{inrg_stage}}{a logical vector}
    \item{\code{enneking_msts_metastasis}}{a logical vector}
    \item{\code{cog_liver_stage}}{a logical vector}
    \item{\code{child_pugh_classification}}{a logical vector}
    \item{\code{metastasis_at_diagnosis_site}}{a logical vector}
    \item{\code{cog_rhabdomyosarcoma_risk_group}}{a logical vector}
    \item{\code{primary_gleason_grade}}{a character vector}
    \item{\code{inpc_grade}}{a logical vector}
    \item{\code{irs_group}}{a logical vector}
    \item{\code{medulloblastoma_molecular_classification}}{a logical vector}
    \item{\code{wilms_tumor_histologic_subtype}}{a logical vector}
    \item{\code{residual_disease}}{a logical vector}
    \item{\code{weiss_assessment_score}}{a logical vector}
    \item{\code{tumor_focality}}{a logical vector}
    \item{\code{ann_arbor_b_symptoms}}{a character vector}
    \item{\code{figo_stage}}{a character vector}
    \item{\code{burkitt_lymphoma_clinical_variant}}{a logical vector}
    \item{\code{days_to_best_overall_response}}{a logical vector}
    \item{\code{inss_stage}}{a logical vector}
    \item{\code{supratentorial_localization}}{a logical vector}
    \item{\code{ishak_fibrosis_score}}{a logical vector}
    \item{\code{tumor_confined_to_organ_of_origin}}{a logical vector}
    \item{\code{gleason_grade_group}}{a logical vector}
    \item{\code{laterality}}{a logical vector}
    \item{\code{ajcc_clinical_stage}}{a character vector}
    \item{\code{cog_neuroblastoma_risk_group}}{a logical vector}
    \item{\code{metastasis_at_diagnosis}}{a logical vector}
    \item{\code{enneking_msts_tumor_site}}{a logical vector}
    \item{\code{secondary_gleason_grade}}{a character vector}
    \item{\code{best_overall_response}}{a logical vector}
    \item{\code{ann_arbor_pathologic_stage}}{a logical vector}
    \item{\code{ann_arbor_extranodal_involvement}}{a character vector}
    \item{\code{method_of_diagnosis}}{a logical vector}
    \item{\code{mitosis_karyorrhexis_index}}{a logical vector}
    \item{\code{ajcc_clinical_n}}{a character vector}
    \item{\code{ajcc_clinical_t}}{a character vector}
    \item{\code{inpc_histologic_group}}{a logical vector}
    \item{\code{masaoka_stage}}{a character vector}
    \item{\code{micropapillary_features}}{a logical vector}
    \item{\code{igcccg_stage}}{a logical vector}
    \item{\code{tumor_regression_grade}}{a logical vector}
    \item{\code{first_symptom_prior_to_diagnosis}}{a logical vector}
    \item{\code{enneking_msts_grade}}{a logical vector}
    \item{\code{gastric_esophageal_junction_involvement}}{a logical vector}
    \item{\code{alcohol_days_per_week}}{a logical vector}
    \item{\code{type_of_smoke_exposure}}{a logical vector}
    \item{\code{smoking_frequency}}{a logical vector}
    \item{\code{type_of_tobacco_used}}{a logical vector}
    \item{\code{alcohol_drinks_per_day}}{a logical vector}
    \item{\code{environmental_tobacco_smoke_exposure}}{a logical vector}
    \item{\code{radon_exposure}}{a logical vector}
    \item{\code{pack_years_smoked}}{a numeric vector}
    \item{\code{asbestos_exposure}}{a logical vector}
    \item{\code{tobacco_smoking_quit_year}}{a logical vector}
    \item{\code{tobacco_smoking_status}}{a logical vector}
    \item{\code{tobacco_smoking_onset_year}}{a logical vector}
    \item{\code{respirable_crystalline_silica_exposure}}{a logical vector}
    \item{\code{time_between_waking_and_first_smoke}}{a logical vector}
    \item{\code{coal_dust_exposure}}{a logical vector}
    \item{\code{cause_of_death}}{a logical vector}
    \item{\code{premature_at_birth}}{a logical vector}
    \item{\code{weeks_gestation_at_birth}}{a logical vector}
    \item{\code{project}}{a character vector}
    \item{\code{bcr_patient_barcode}}{a character vector}
    \item{\code{disease}}{a character vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
GDC v32
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
}
\keyword{datasets}

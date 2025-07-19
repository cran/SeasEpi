\docType{data}
\name{Datasets}
\alias{adjacency_matrix}
\alias{data}
\title{Hypothetical Datasets}
\description{
The main function, \code{SeasEpi_Par_Est}, applies the spatiotemporal individual-level model of seasonal infectious disease transmission within the SEIRS framework to a hypothetical dataset. It is compatible with any dataset that follows the required format, which includes two dataframes: \code{data} and \code{adjacency_matrix}, along with relevant parameter inputs. To demonstrate the expected input structure and the function’s practical use, we provide two hypothetical examples of \code{data} and \code{adjacency_matrix}.
}
\section{\code{data}}{
A data frame with 100 rows and 11 columns.

This sample dataset illustrates the required structure for the dataframe used with this package. While the number of rows can vary, each row must represent a single infected individual, and the column names and order must follow the specified format. The example includes individual-level attributes (e.g., age, infection status) as well as area-level information (e.g., socioeconomic status) for 100 individuals, each linked to a postal code. This dataset will serve as input in the example demonstrating the \code{SeasEpi_Par_Est} function.
\describe{
\item{Ave_Postal_Pop}{Average population of each postal code}
\item{AverageAge}{Average age of individuals within each postal code (individual-level data)}
\item{InfectedTime}{Time of infection for each individual, represented as a numerical value from 1 to the end of the pandemic period}
\item{LAT}{Latitude of the postal code}
\item{LONG}{Longitude of the postal code}
\item{Label_NC_shape}{The region number that the postal code belongs to, here assuming the study area is divided into five subregions}
\item{MaleRate}{Rate of males in the population of the postal code (individual-level data)}
\item{NInfected}{Number of infected individuals in the postal code}
\item{SES}{Socioeconomic status indicator of the region to which the postal code belongs (area-level data)}
\item{STI}{Sexually transmitted infection rate of the region that the postal code belongs to (area-level data)}
\item{SymptomRate}{Rate of disease symptoms in the postal code (individual-level data), indicating whether individuals are symptomatic or asymptomatic}
}
}



\section{\code{adjacency_matrix}}{
A 5x5 matrix.

This hypothetical adjacency matrix is provided to illustrate the structure required for use with this package. The matrix used with the package should follow a similar format, maintaining the same layout but allowing for any number of regions. The adjacency matrix  defines the neighborhood relationships between subregions in a hypothetical study area. In this example, it represents a spatial structure with five subregions, where each cell indicates the presence or absence of a connection between the corresponding subregions. The example for the \code{SeasEpi_Par_Est} function will use this matrix as input.
\describe{
\item{V1}{Subregion 1: Represents the first subregion in the region under study}
\item{V2}{Subregion 2: Represents the second subregion in the region under study}
\item{V3}{Subregion 3: Represents the third subregion in the region under study}
\item{V4}{Subregion 4: Represents the fourth subregion in the region under study}
\item{V5}{Subregion 5: Represents the fifth subregion in the region under study}
\item{Value}{Each cell in the matrix (e.g., between subregion 1 and subregion 2) represents the connection (typically 0 or 1) between the two subregions, where 1 indicates they are neighbors and 0 indicates they are not.}
}
}



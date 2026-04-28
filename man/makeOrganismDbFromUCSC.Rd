\name{makeOrganismDbFromUCSC}
\alias{makeOrganismDbFromUCSC}

\title{
  Make a OrganismDb object from annotations available at the
  UCSC Genome Browser
}
\description{
  The \code{makeOrganismDbFromUCSC} function allows the user
  to make a \link{OrganismDb} object from transcript annotations
  available at the UCSC Genome Browser.
}
\usage{
makeOrganismDbFromUCSC(
        genome="hg19",
        tablename="knownGene",
        transcript_ids=NULL,
        circ_seqs=NULL,
        goldenPath.url=getOption("UCSC.goldenPath.url"))
}
\arguments{
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{\link[rtracklayer]{ucscGenomes}()[ , "db"]}.
    For example: \code{"hg19"}.}
  \item{tablename}{name of the UCSC table containing the transcript
    annotations to retrieve. Use the \code{supportedUCSCtables} utility
    function to get the list of supported tables. Note that not all
    tables are available for all genomes.}
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{OrganismDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{goldenPath.url}{use to specify the location of an
    alternate UCSC Genome Browser.}
}
\details{
  \code{makeOrganismDbFromUCSC} is a convenience function that feeds
  data from the UCSC source to the lower level \code{\link{OrganismDb}}
  function.
  See \code{?\link{makeOrganismDbFromBiomart}} for a similar function
  that feeds data from a BioMart database.
}

\value{A \link{OrganismDb} object.}

\author{
  M. Carlson
}

\seealso{
  \itemize{
    \item \code{\link{makeOrganismDbFromBiomart}} for convenient ways to make a
          \link{OrganismDb} object from BioMart online resources.

    \item \code{\link[rtracklayer]{ucscGenomes}} in the \pkg{rtracklayer}
          package.

    \item The \link{OrganismDb} class.
  }
}

\examples{
\dontrun{
## Display the list of genomes available at UCSC:
library(rtracklayer)
library(RMariaDB)
ucscGenomes()[ , "db"]

## Display the list of tables supported by makeOrganismDbFromUCSC():
supportedUCSCtables()

\dontrun{
## Retrieving a full transcript dataset for Yeast from UCSC:
odb1 <- makeOrganismDbFromUCSC(genome="sacCer2", tablename="ensGene")
}

## Retrieving an incomplete transcript dataset for Mouse from UCSC
## (only transcripts linked to Entrez Gene ID 22290):
transcript_ids <- c(
    "uc009uzf.1",
    "uc009uzg.1",
    "uc009uzh.1",
    "uc009uzi.1",
    "uc009uzj.1"
)

odb2 <- makeOrganismDbFromUCSC(genome="mm9", tablename="knownGene",
                          transcript_ids=transcript_ids)
odb2
}
}

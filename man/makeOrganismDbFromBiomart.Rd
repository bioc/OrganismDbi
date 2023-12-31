\name{makeOrganismDbFromBiomart}
\alias{makeOrganismDbFromBiomart}

\title{
  Make a OrganismDb object from annotations available on a
  BioMart database
}
\description{
  The \code{makeOrganismDbFromBiomart} function allows the user
  to make a \link{OrganismDb} object from transcript annotations
  available on a BioMart database.  This object has all the benefits of
  a TxDb, plus an associated OrgDb and GODb object.
}
\usage{
makeOrganismDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                    dataset="hsapiens_gene_ensembl",
                    transcript_ids=NULL,
                    circ_seqs=NULL,
                    filter="",
                    id_prefix="ensembl_",
                    host="https://www.ensembl.org",
                    port,
                    miRBaseBuild=NA,
                    keytype = "ENSEMBL",
                    orgdb = NA)

}
\arguments{
  \item{biomart}{which BioMart database to use.
    Get the list of all available BioMart databases with the
    \code{\link[biomaRt]{listMarts}} function from the biomaRt
    package. See the details section below for a list of BioMart
    databases with compatible transcript annotations.}
  \item{dataset}{which dataset from BioMart. For example:
    \code{"hsapiens_gene_ensembl"}, \code{"mmusculus_gene_ensembl"},
    \code{"dmelanogaster_gene_ensembl"}, \code{"celegans_gene_ensembl"},
    \code{"scerevisiae_gene_ensembl"}, etc in the ensembl database.
    See the examples section below for how to discover which datasets
    are available in a given BioMart database.}
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{TxDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.  This \link{TxDb} object
    will be embedded in the resulting \link{OrganismDb} object.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{filter}{Additional filters to use in the BioMart query. Must be
    a named list. An example is \code{filter=as.list(c(source="entrez"))}}
  \item{host}{The host URL of the BioMart. Defaults to www.ensembl.org.}
  \item{port}{Deprecated: The port to use in the HTTP communication
    with the host.}
  \item{id_prefix}{Specifies the prefix used in BioMart attributes. For
    example, some BioMarts may have an attribute specified as
    \code{"ensembl_transcript_id"} whereas others have the same attribute
    specified as \code{"transcript_id"}. Defaults to \code{"ensembl_"}.}
  \item{miRBaseBuild}{specify the string for the appropriate build
    Information from mirbase.db to use for microRNAs.  This can be
    learned by calling \code{supportedMiRBaseBuildValues}.  By default,
    this value will be set to \code{NA}, which will inactivate the
    \code{microRNAs} accessor.}
  \item{keytype}{This indicates the kind of key that this database will
    use as a foreign key between it's TxDb object and it's OrgDb
    object. So basically whatever the column name is for the foreign key
    from your OrgDb that your TxDb will need to map it's GENEID on to.
    By default it is "ENSEMBL" since the GENEID's for most biomaRt based
    TxDbs will be ensembl gene ids and therefore they will need to map
    to ENSEMBL gene mappings from the associated OrgDb object.}
  \item{orgdb}{By default, \code{makeOrganismDbFromBiomart} will use the
    taxonomyID from your txdb to lookup an appropriate matching
    \code{OrgDb} object but using this you can supply a different
    \code{OrgDb} object.}
}

\details{
  \code{makeOrganismDbFromBiomart} is a convenience function that feeds
  data from a BioMart database to the lower level
  \code{\link{OrganismDb}} constructor.
  See \code{?\link{makeOrganismDbFromUCSC}} for a similar function
  that feeds data from the UCSC source.

  The \code{listMarts} function from the \pkg{biomaRt} package can be
  used to list all public BioMart databases.
  Not all databases returned by this function contain datasets that
  are compatible with (i.e. understood by) \code{makeOrganismDbFromBiomart}.
  Here is a list of datasets known to be compatible (updated on Sep 24, 2014):
  \itemize{
    \item All the datasets in the main Ensembl database:
          use \code{biomart="ensembl"}.

    \item All the datasets in the Ensembl Fungi database:
          use \code{biomart="fungi_mart_XX"} where XX is the release
          version of the database e.g. \code{"fungi_mart_22"}.

    \item All the datasets in the Ensembl Metazoa database:
          use \code{biomart="metazoa_mart_XX"} where XX is the release
          version of the database e.g. \code{"metazoa_mart_22"}.

    \item All the datasets in the Ensembl Plants database:
          use \code{biomart="plants_mart_XX"} where XX is the release
          version of the database e.g. \code{"plants_mart_22"}.

    \item All the datasets in the Ensembl Protists database:
          use \code{biomart="protists_mart_XX"} where XX is the release
          version of the database e.g. \code{"protists_mart_22"}.

    \item All the datasets in the Gramene Mart:
          use \code{biomart="ENSEMBL_MART_PLANT"}.
  }
  Not all these datasets have CDS information.
}

\value{A \link{OrganismDb} object.}

\author{
  M. Carlson
}

\seealso{
  \itemize{
    \item \code{\link{makeOrganismDbFromUCSC}} for convenient ways to make a
          \link{OrganismDb} object from UCSC online resources.

    \item The \code{\link[biomaRt]{listMarts}}, \code{\link[biomaRt]{useMart}},
          and \code{\link[biomaRt]{listDatasets}} functions in the
          \pkg{biomaRt} package.

    \item The \code{\link{supportedMiRBaseBuildValues}} function for
          listing all the possible values for the \code{miRBaseBuild}
          argument.

    \item The \link{OrganismDb} class.

  }
}

\examples{
## Discover which datasets are available in the "ensembl" BioMart
## database:
library(biomaRt)
mart <- useEnsembl("ensembl")
datasets <- listDatasets(mart)
head(datasets)

## Retrieving an incomplete transcript dataset for Human from the
## "ensembl" BioMart database:
transcript_ids <- c(
    "ENST00000013894",
    "ENST00000268655",
    "ENST00000313243",
    "ENST00000435657",
    "ENST00000384428",
    "ENST00000478783"
)
odb <- makeOrganismDbFromBiomart(transcript_ids=transcript_ids)
odb  # note that these annotations match the GRCh38 genome assembly

if (interactive()) {
  ## Now what if we want to use another mirror?  We might make use of the
  ## new host argument.  But wait!  If we use biomaRt, we can see that
  ## this host has named the mart differently!
  listMarts(host="https://useast.ensembl.org")

  ## Therefore we must also change the name passed into the "mart"
  ## argument thusly:
  makeOrganismDbFromBiomart(
      biomart="ENSEMBL_MART_ENSEMBL",
      transcript_ids=transcript_ids,
      host="https://useast.ensembl.org"
  )
}

}

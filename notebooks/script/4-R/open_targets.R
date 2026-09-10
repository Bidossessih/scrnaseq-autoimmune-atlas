# Install relevant library for HTTP requests
library(httr)

# Set gene_id variable for AR (androgen receptor)
gene_id <- "ENSG00000169083"

# Build query string to get general information about AR and genetic constraint and tractability assessments 
query_string = "
  query KnownDrugsQuery(
  $efoId: String!
  $cursor: String
  $freeTextQuery: String
  $size: Int = 10
) {
  disease(efoId: 'EFO_0005140') {
    id
    knownDrugs() {
      count
      cursor
      rows {
        phase
        status
        urls {
          name
          url
        }
        disease {
          id
          name
        }
        drug {
          id
          name
          mechanismsOfAction {
            rows {
              actionType
              targets {
                id
              }
            }
          }
        }
        urls {
          url
          name
        }
        drugType
        mechanismOfAction
        target {
          id
          approvedName
          approvedSymbol
        }
      }
    }
  }
}

"

# Set base URL of GraphQL API endpoint
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"

# Set variables object of arguments to be passed to endpoint
variables <- list("ensemblId" = gene_id)

# Construct POST request body object with query string and variables
post_body <- list(query = query_string)

# Perform POST request
r <- POST(url=base_url, body=post_body, encode='json')

# Print data to RStudio console
print(content(r)$data)

# Response to Reviewers

We thank the reviewers for their thoughtful comments that have helped improve our manuscript. Below, we address each comment and detail the changes made to the manuscript.

## Changes Already Implemented

1. **Fixed citations**:
   - Citation [10] (Tarca2013) has been appropriately moved to the KEGGdzPathwaysGEO package description
   - Reference [15] (Yosef2011) has been added to the ANAT tool citation

2. **Clarified the meaning of "rank" in the evaluation**:
   - Added explanation in line 91 that "rank refers to the position of the true associated pathway in the list of all pathways sorted by their enrichment significance, with lower ranks indicating better performance (i.e., rank 1 means the true pathway was identified as the most significant pathway)"

3. **Web Server Reliability**:
   - Fixed web server issues by implementing dedicated download endpoints with proper headers
   - Added request logging and health check endpoint for improved monitoring
   - Updated Docker configuration to ensure consistent server operation
   - The web server is now continuously available at https://peanut.cs.tau.ac.il/

4. **Code Repository and DOI**:
   - Created GitHub repository with complete code and documentation
   - Archived code on Zenodo and obtained a DOI for long-term preservation

## Changes Still Needed

1. **Add detailed comparison between PEANUT and NGSEA methodologies**:
   - Need to add a new paragraph titled "Comparison between NGSEA and PEANUT" in the Methods section
   - Explain how PEANUT's global network propagation approach differs from NGSEA's localized neighbor-based approach
   - Discuss how these methodological differences affect performance

2. **Explain choice of permutation tests**:
   - Add justification for using 10,000 permutations as a balance between statistical power and computational efficiency

3. **Add Code Availability Section**:
   - Add a new "Code and Data Availability" section to the paper
   - Include GitHub repository URL and Zenodo DOI

## Reviewer 1 Comments

The reviewer noted that our "experimental results demonstrate the advantages of a network propagation based approach over existing gene expression analysis methods." We appreciate this positive assessment of our work.

## Reviewer 2 Comments

### Comment 1: Comparison with NGSEA
> "Although the authors show in Fig. 1 that PEANUT produced better results than GSEA, ABS GSEA and NGSEA, the NGSEA paper showed the opposite results regarding the comparison between GSEA and NGSEA... it would be better to discuss their differences in more detail."

**Response**: We will add a dedicated paragraph titled "Comparison between NGSEA and PEANUT" that explains in detail the methodological differences between PEANUT and NGSEA. We will discuss how PEANUT uses global network propagation through Random Walk with Restart, which diffuses information throughout the entire network, while NGSEA only considers immediate neighbors. We will also address the potential reasons for the different results between our implementation and the original NGSEA paper.

### Comment 2: Web Server Availability
> "Unfortunately, there were no responses from the web server when I accessed it, so I could not evaluate it for its user-friendliness, etc."

**Response**: We apologize for the server issues the reviewer experienced. We have thoroughly fixed the web server issues by implementing dedicated download endpoints with proper headers, adding request logging and a health check endpoint for improved monitoring, and updating the Docker configuration to ensure consistent operation. The web server is now continuously available at https://peanut.cs.tau.ac.il/.

### Comment 3: Define Meaning of Rank
> "It would be better to add the meaning of the rank to easily understand the evaluation results."

**Response**: We have added clarification about what "rank" means in our evaluation context in the Results section, explaining that it refers to the position of the true associated pathway in the sorted list, with lower ranks indicating better performance.

### Comment 4: Citation Corrections
> "The citation [10] would be better to move to the sentence describing KGGdzPathwayGEO. The reference 15 should be added to the ANAT tool paragraph."

**Response**: We have made both of these changes as suggested.

### Comment 5: Permutation Tests Justification
> "It would be nice to explain how the number of the iterations of 10,000 for the permutation tests is determined regarding the effectiveness and the computational resource."

**Response**: We will add an explanation regarding our choice of 10,000 permutations, noting that it represents a balance between achieving reliable statistical significance and maintaining reasonable computational efficiency.

## Additional Changes Required by the Journal

The journal requires us to archive our code on a site like Zenodo to obtain a DOI. We have completed this requirement and will add a "Code and Data Availability" section to the paper containing the repository URL and DOI information. 
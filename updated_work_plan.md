# Updated Work Plan for PEANUT Paper Revision

## Immediate Actions (Next 1-2 Days)

1. **Add NGSEA Comparison** ⚠️ HIGH PRIORITY
   - Add a new paragraph titled "Comparison between NGSEA and PEANUT" to the Methods section
   - Explain how PEANUT uses global network propagation while NGSEA only considers immediate neighbors
   - Address the different results between our implementation and the original NGSEA paper
   - Timeline: Complete within 1-2 days

2. **Explain Permutation Test Choice** ⚠️ HIGH PRIORITY
   - Add justification for using 10,000 permutations 
   - Explain the balance between statistical power and computational efficiency
   - Run multiple steps comparing
   - Timeline: Complete within 1 day

3. **Add Code Availability Section**
   - Add a new "Code and Data Availability" section to the paper
   - Include GitHub repository URL and Zenodo DOI
   - Timeline: Complete within 2 days

## Final Paper Revision (Within 4 Days)

4. **Additional Paper Edits** 
   - Conduct final proofreading of all changes made
   - Check for any remaining issues or inconsistencies
   - Ensure all reviewer comments have been fully addressed
   - Timeline: Complete within 4 days

5. **Supplementary Materials**
   - Prepare clean supplementary materials in PDF/Word format
   - Ensure all files are properly formatted according to journal requirements
   - Timeline: Complete within 4 days

## Submission Preparation (Within 6 Days)

6. **Final Response Document**
   - Update the revision response document with:
     - Confirmation of all implemented changes
     - Web server functionality status
     - Code repository and DOI information
   - Timeline: Complete within 6 days

7. **Final Submission Package**
   - Prepare all files for submission:
     - Revised manuscript with changes highlighted
     - Response to reviewers document
     - All required LaTeX files
     - Supplementary materials
   - Timeline: Submit within 6 days

## Post-submission (Ongoing)

8. **Monitoring and Maintenance**
   - Continue monitoring web server for any issues
   - Address any questions from editors or reviewers promptly
   - Timeline: Throughout the review process

## Completed Tasks ✓

- Fixed citation issues:
  - Moved citation [10] (Tarca2013) to the KEGGdzPathwaysGEO package description
  - Added reference [15] (Yosef2011) to the ANAT tool citation
- Added explanation of "rank" in the Results section
- Created initial revision response document
- Fixed web server reliability issues:
  - Added dedicated download endpoints with proper headers
  - Implemented request logging and health check endpoint
  - Updated Dockerfile with necessary changes
- Created GitHub repository and obtained Zenodo DOI 
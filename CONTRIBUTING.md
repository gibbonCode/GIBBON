# Contribution Guidelines

## Reporting issues

- **Search for existing issues.** Please check to see if someone else has reported the same issue.
- **Share as much information as possible.** Include operating system and version, browser and version. Also, include steps to reproduce the bug.

## Project Setup
Refer to the [README](README.md).

## Code Style

### Variable Naming
Not all current code follows the conventions below but these will be followed for future developments. 
- `lowerCamelCase` General variables
- `UpperCamelCase` Functions
- Maximize the use  of semantic and descriptive variables names (e.g. `faceIndices` not `fcInd` or `fi`). Avoid abbreviations except in cases of industry wide usage. In some cases non-descriptive and short variable names are exceptable for instance vertices (points), faces, edges, colors and logic arrays may be denoted `V`, `F`, `E`, `C`, `L`. Furthermore, if a mathematrical symbol or letter is commonly used for some entity it may be acceptable to use short names e.g. coordinates may be referred to as `X`, `Y` and `Z` and image coordinates of indices may be referred to as `I`, `J` and `K`. In some cases the use of capital or non-capital letters refers to tensors/matrices/arrays/sets and scalars/components/subsets respectively, e.g. a multitude of scalars `c` may be contained within an array or matrix `C`, or a cell array `D` may contain individual entries referred to as `d`. 

## Testing
For the moment the DEMO_ and HELP_ files may serves as a test suite. 

## Pull requests
- Try not to pollute your pull request with unintended changes – keep them simple and small. If possible, squash your commits.
- Try to share how your code has been tested before submitting a pull request.
- If your PR resolves an issue, include **closes #ISSUE_NUMBER** in your commit message (or a [synonym](https://help.github.com/articles/closing-issues-via-commit-messages)).
- Review
    - If your PR is ready for review, another contributor will be assigned to review your PR
    - The reviewer will accept or comment on the PR. 
    - If needed address the comments left by the reviewer. Once you're ready to continue the review, ping the reviewer in a comment.
    - Once accepted your code will be merged to `master`

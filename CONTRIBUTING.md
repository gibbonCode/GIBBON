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

## Documentation
The documentation is found in the `GIBBON/docs` folder in the form of `HELP_functionName.m` files and `DEMO_demoName.m` files. To create documentation for a particular function e.g. named `functionName.m` one can add a `HELP_functionName.m` file containing the documentation. Please use [MATLAB's "publishing markup"](https://uk.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html) (it may be useful to copy an existing help/demo file to get started), which will allow the documentation files to be converted to HTML files to serve as online documentation. Once documentation is ready to be added to the website https://www.gibboncode.org/Documentation/ you can run:
```matlab
gpublish HELP_functionName
```
This will "publish" the documentation in the form of HTML files in the `GIBBON/docs/html` folder. The content there will then be automatically added to, and rendered on, the website when the website is updated. 

## Testing
For the moment the DEMO_ and HELP_ files may serves as a test suite. The `testGibbon` function can be used to run (and "publish" HTML if needed) these tests/demos automatically. 

## Pull requests
- Try not to pollute your pull request with unintended changes – keep them simple and small. If possible, squash your commits.
- Try to share how your code has been tested before submitting a pull request.
- If your PR resolves an issue, include **closes #ISSUE_NUMBER** in your commit message (or a [synonym](https://help.github.com/articles/closing-issues-via-commit-messages)).
- Review
    - If your PR is ready for review, another contributor will be assigned to review your PR
    - The reviewer will accept or comment on the PR. 
    - If needed address the comments left by the reviewer. Once you're ready to continue the review, ping the reviewer in a comment.
    - Once accepted your code will be merged to `master`

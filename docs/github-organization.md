# GitHub Best Practices Guide

## Branching Strategy

### Branch Naming
- Use descriptive, kebab-case names that summarize the purpose
- Prefix with feature/, fix/, docs/, test/, refactor/, or chore/
- Examples:
  - `feature/add-user-authentication`
  - `fix/issue-42-login-bug`
  - `docs/update-readme`

### Branch Management
- Create a new branch for each feature or bugfix
- Keep branches focused on a single purpose
- Keep branches up-to-date with main/master regularly
- Delete branches after they've been merged

## Pull Requests (PRs)

### Creating PRs
- Open a PR for any changes to the main/master branch
- Keep PRs small and focused on a single feature/fix
- Include a clear, descriptive title
- Use the PR template if available
- Link to any related issues

### PR Description
- Describe the purpose of the changes
- Document any breaking changes
- Include testing instructions
- Add screenshots/GIFs for UI changes
- Note any dependencies or special setup required

### Code Review
- Request reviews from relevant team members
- Use GitHub's "Request Changes" feature for required changes
- Resolve all comments before merging
- Use GitHub's "Suggestion" feature for small code suggestions
- Include Copilot review when appropriate

## Data Management

### Large Files
- **GitHub File Size Limit**: Keep individual files under 100MB
- **For files < 700MB**:
  - Consider using Git LFS (Large File Storage)
  - Document setup instructions in README
- **For files > 700MB**:
  - Store on Google Drive or similar service
  - Include clear download instructions in README
  - Specify exact directory structure needed
  - Consider including a subsampled version in the repo for testing

### Data Organization
- Keep data separate from code
- Document data sources and processing steps
- Include a data dictionary for complex datasets
- Use `.gitignore` to exclude large or sensitive files

### Best Practices
- Never commit sensitive information (API keys, credentials)
- Use environment variables for configuration
- Include a `.env.example` file with required environment variables
- Document data dependencies and setup steps in README.md


# prepub_supplementary
Paper supplements before publication

# Linking to public version

From [StackOverflow](https://stackoverflow.com/questions/7983204/having-a-private-branch-of-a-public-repo-on-github)

1. Duplicate your repo.
2. Make the duplicated repo a private one on GitHub.
3. Clone the private repo to your machine
4. Add a remote to your public repo (git remote add public git@github.com:...)
5. Push branches with commits intended for your public repo to that new public remote. (make sure you don't accidentally commit private-only code)
6. You can bring in changes to your public repo using 'git fetch public' and then merge them locally and push to your private repo (origin remote).
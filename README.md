# master-project-backup

## back up procedure
```bash=
sh update.sh
git add *
git commit
git push
```

## version when in branch
```bash=
git checkout branch
git add *
git commit
git push
git checkout master
git merge branch
git push origin master (not sure)
```
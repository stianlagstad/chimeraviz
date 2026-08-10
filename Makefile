GIT := \git

.PHONY: sync-from-bioc publish-to-bioc

sync-from-bioc:
	$(GIT) fetch upstream
	$(GIT) checkout master
	$(GIT) merge --ff-only upstream/devel
	$(GIT) push origin master

publish-to-bioc:
	$(GIT) checkout master
	$(GIT) push origin master
	$(GIT) push upstream master:devel

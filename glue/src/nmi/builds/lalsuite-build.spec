# Metronome submit file for CBC builds

# metadata for the Metronome DB
description = build of lalsuite from repo=$(GIT_REPO), branch=$(GIT_BRANCH), id=$(GIT_ID)
project = cbc
project_version = $(GIT_BRANCH)
component = lalsuite
component_version = $(GIT_ID)
run_type = build

# build parameters
platforms = x86_64_cent_5.3
inputs = $(HARNESS_INPUT_SPEC_FILE)
remote_task = remote_task.sh
remote_post = remote_post.py
remote_post_args = --verbose

# send the submitting user an email if any task fails
notify = $(USER)@syr.edu
notify_fail_only = true

# stream "live" stdout/err from tasks back to the submit node
stream_user_io = true

# we set these attributes so that they get propagated to the
# environment of our input spec files and/or remote_* scripts above
git_repo = $(GIT_REPO)
git_id = $(GIT_ID)
git_branch = $(GIT_BRANCH)

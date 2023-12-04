dircreater <- function(directory, recursive=TRUE) {
  if (!file.exists(directory)){
    dir.create(directory, recursive= recursive)
  }
  else {print('Directory does already exist')}
} 
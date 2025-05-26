function normalizedImages = normalizeArray(images)
normalizedImages = (images - min(images(:))) / (max(images(:)) - min(images(:)));
function val = isBlockDimsWellFormated(blockDim)

warning("deprecated:isBlockDimsWellFormated",...
    "isBlockDimsWellFormated is deprecated, use isBlockDimsWellFormatted instead.")

val = isBlockDimsWellFormatted(blockDim);
end
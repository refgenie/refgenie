# TO be imported from refget package when it is finished

def parse_fasta(self, fa_file):
    _LOGGER.info("Hashing {}".format(fa_file))
    try:
        fa_object = pyfaidx.Fasta(fa_file)
    except pyfaidx.UnsupportedCompressionFormat:
        # pyfaidx can handle bgzip but not gzip; so we just hack it here and
        # unzip the file for checksumming, then rezip it for the rest of the
        # asset build.
        # TODO: streamline this to avoid repeated compress/decompress
        os.system("gunzip {}".format(fa_file))
        fa_file_unzipped = fa_file.replace(".gz", "")
        fa_object = pyfaidx.Fasta(fa_file_unzipped)
        os.system("gzip {}".format(fa_file_unzipped))
    return fa_object

def fasta_checksum(self, fa_file):
    """
    Just calculate checksum of fasta file without loading it.
    """
    fa_object = parse_fasta(fa_file)
    content_checksums = {}
    for k in fa_object.keys():
        content_checksums[k] = self.checksum_function(str(fa_object[k]))
    collection_string = ";".join([":".join(i) for i in content_checksums.items()])
    collection_checksum = self.load_seq(collection_string)
    return collection_checksum, content_checksums
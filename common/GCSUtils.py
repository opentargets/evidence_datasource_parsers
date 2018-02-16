from settings import Config
import logging
import google.cloud.storage as gcs
import google.cloud.exceptions as gce
import tempfile
import os


class GCSBucketManager(object):

    def __init__(self):
        print("init GCSBucketManager")
        super(GCSBucketManager, self).__init__()
        self._logger = logging.getLogger(__name__)
        self.gcs_client = gcs.Client(Config.GOOGLE_DEFAULT_PROJECT)
        self.bucket = None
        self._logger.info("This is done")

        try:
            for bucket in self.gcs_client.list_buckets():
                print(bucket)
            self.bucket = self.gcs_client.get_bucket(Config.GOOGLE_BUCKET_EVIDENCE_INPUT)
        except gce.NotFound:
            self._logger.error('Please set export GOOGLE_APPLICATION_CREDENTIALS=/path/to/key.json')
            print('Sorry, that bucket does not exist!')


    def get_gcs_filename(self, filename):

        blob = self.bucket.get_blob(filename)
        file_obj = tempfile.NamedTemporaryFile(delete=False)
        name = file_obj.name
        blob.download_to_file(file_obj)
        return name

    def upload_blob_from_string(self, data, filename):

        blob = self.bucket.blob(filename)
        blob.upload_from_string(data, content_type='text/plain')

    def download_blob_as_string(self, filename):

        print("download %s"%filename)
        blob = self.bucket.get_blob(filename)
        print(blob.exists())
        raw = blob.download_as_string()
        print(raw)
        return raw

    def upload_gcs_file(self, filename):

        tmpdir = tempfile.mkdtemp()

        # Ensure the file is read/write by the creator only
        saved_umask = os.umask(0o077)
        path = os.path.join(tmpdir, filename)

        try:
            with open(path, "w") as tmp:
                tmp.write("secrets!")
        except IOError as e:
            print
            'IOError'
        else:
            os.remove(path)
        finally:
            os.umask(saved_umask)
            os.rmdir(tmpdir)


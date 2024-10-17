import boto3
import os
import tarfile
from botocore.exceptions import ClientError        
from multiprocessing.pool import ThreadPool

import cht_utils.fileops as fo

class S3Session:
    # Helper class for cloud functions

    def __init__(self, access_key, secret_key, region):  

        self.ready = True
        
        # Create a session using your AWS credentials (or configure it in other ways)
        session = boto3.Session(
            aws_access_key_id=access_key,
            aws_secret_access_key=secret_key,
            region_name=region
        )
        # Create an S3 client
        self.s3_client = session.client('s3')

    def upload_file(self, bucket_name, file, s3_folder, quiet=True):
        s3_key = os.path.join(s3_folder, os.path.basename(file)).replace('\\', '/')
        self.s3_client.upload_file(file, bucket_name, s3_key)
        if not quiet:
            print("Uploaded " + os.path.basename(file))

    def download_file(self, bucket_name, s3_folder, file, local_folder, quiet=True):
        s3_key = os.path.join(s3_folder, os.path.basename(file)).replace('\\', '/')
        local_path = os.path.join(local_folder, os.path.basename(file))
        self.s3_client.download_file(bucket_name, s3_key, file, local_path)
        if not quiet:
            print("Downloaded " + os.path.basename(file))

    def download_files(self, bucket_name, key_list, file_list, quiet=True):
        # key_list is a list of keys in the bucket
        # file_list is a list of local file paths
        pool = ThreadPool()
        pool.starmap(self.download_file_parallel, [(bucket_name, key, file, quiet) for key, file in zip(key_list, file_list)])

    def download_file_parallel(self, bucket_name, s3_key, file, quiet=True):
        self.s3_client.download_file(bucket_name, s3_key, file)
        if not quiet:
            print("Downloaded " + os.path.basename(file))

    def delete_file(self, bucket_name, s3_folder, file, quiet=True):
        s3_key = os.path.join(s3_folder, os.path.basename(file)).replace('\\', '/')
        self.s3_client.delete_object(Bucket=bucket_name, Key=s3_key)
        if not quiet:
            print("Deleted " + os.path.basename(file))

    def make_folder(self, bucket_name, s3_folder, quiet=True):
        self.s3_client.put_objects(bBucket=bucket_name, Key=(s3_folder + '/'))
        if not quiet:
            print("Made folder: " + s3_folder)

    def upload_folder(self, bucket_name, local_folder, s3_folder, parallel=True, quiet=True):
        local_folder = local_folder.replace('\\\\','\\')
        local_folder = local_folder.replace('\\','/')
        # Recursively list all files
        flist = list_all_files(local_folder)
        if parallel:
            pool = ThreadPool()
            pool.starmap(upf, [(file, local_folder, s3_folder, bucket_name, self.s3_client, quiet) for file in flist])
        else:
            for file in flist:
                upf(file, local_folder, s3_folder, bucket_name, self.s3_client, quiet)

    def download_folder(self, bucket_name, s3_folder, local_folder, parallel=True, quiet=True):
        fo.mkdir(local_folder)
        objects = self.s3_client.list_objects(Bucket=bucket_name, Prefix=s3_folder)
        if "Contents" in objects:
            # if parallel:
            #     # There may be a faster way for this with imap and global variables?
            #     pool = ThreadPool()
            #     pool.starmap(self.download_file, [(bucket_name, s3_folder, object['Key'], local_folder, quiet) for object in objects['Contents']])
            # else:
            for object in objects['Contents']:
                s3_key = object['Key']
                local_path = os.path.join(local_folder, os.path.basename(s3_key))
                self.s3_client.download_file(bucket_name, s3_key, local_path)
                if not quiet:
                    print("Downloaded " + os.path.basename(s3_key))

    def delete_folder(self, bucket_name, s3_folder):
        if s3_folder[-1] != "/":
             s3_folder = s3_folder + "/"
        objects = self.s3_client.list_objects(Bucket=bucket_name, Prefix=s3_folder)
        if "Contents" in objects:
            for object in objects['Contents']:
                self.s3_client.delete_object(Bucket=bucket_name, Key=object['Key'])

    def list_folders(self, bucket_name, s3_folder):
        if s3_folder[-1] != "/":
             s3_folder = s3_folder + "/"
        folders = []
        paginator = self.s3_client.get_paginator('list_objects_v2')
        iterator = paginator.paginate(Bucket=bucket_name, Prefix=s3_folder, Delimiter='/')
        for page in iterator:
            for subfolder in page.get('CommonPrefixes', []):
                subfolder_name = subfolder['Prefix'].rstrip('/').split('/')[-1]
                folders.append(subfolder_name)
        return folders 

    def list_files(self, bucket_name, s3_folder):
        paginator = self.s3_client.get_paginator('list_objects_v2')
        
        all_files = []
        
        for page in paginator.paginate(Bucket=bucket_name, Prefix=s3_folder):
            if 'Contents' in page:
                for obj in page['Contents']:
                    all_files.append(obj['Key'])
        
        return all_files

    def download_and_extract_tgz(self, bucket_name, s3_folder, local_folder):
        """
        Download and extract a .tgz file from S3.
        """
        local_tgz_path = os.path.join('/tmp', os.path.basename(s3_folder))
        
        # Download the .tgz file
        self.s3_client.download_file(bucket_name, s3_folder, local_tgz_path)
        
        # Extract the .tgz file
        with tarfile.open(local_tgz_path, "r:gz") as tar:
            tar.extractall(path=local_folder)
        
        # Clean up the downloaded .tgz file
        os.remove(local_tgz_path)

    def check_file_exists(self, bucket_name, s3_key):
        try:
            self.s3_client.head_object(Bucket=bucket_name, Key=s3_key)
            return True
        except ClientError as e:
            if e.response['Error']['Code'] == '404':
                return False
            else:
                raise

    def check_folder_exists(self, bucket_name, s3_key):
        response = self.s3_client.list_objects_v2(Bucket=bucket_name, Prefix=s3_key, Delimiter='/')
        # Check if any items are returned
        if 'CommonPrefixes' in response:
            return True
        else:
            return False

def list_all_files(src):
    # Recursively list all files and folders in a folder
    import pathlib
    pth = pathlib.Path(src)
    pthlst = list(pth.rglob("*"))
    lst = []
    for f in pthlst:
        if f.is_file():
            lst.append(str(f))
    return lst        

def upf(file, local_folder, s3_folder, bucket_name, s3_client, quiet):
    file1 = file.replace('\\','/')
    file1 = file1.replace(local_folder,'')
    s3_key = s3_folder + file1
    s3_client.upload_file(file, bucket_name, s3_key)
    if not quiet:
        print("Uploaded " + file)

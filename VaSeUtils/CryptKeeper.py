from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.asymmetric import rsa, padding
from cryptography.hazmat.primitives import serialization, hashes
from cryptography.fernet import Fernet


class Locker:
    def __init__(self):
        self.private_key = ""
        self.public_key = ""
        self.encrypted_fernet = ""
        return

    def make_key_pair(self):
        self.private_key = rsa.generate_private_key(
                public_exponent=65537,
                key_size=4096,
                backend=default_backend()
                )
        self.public_key = self.private_key.public_key()
        return

    def store_key_pair(self, key_out):
        pem = self.private_key.private_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PrivateFormat.PKCS8,
                encryption_algorithm=serialization.NoEncryption()
                )

        pem2 = self.public_key.public_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PublicFormat.SubjectPublicKeyInfo
                )

        with open(f"{key_out}.pvk", "wb") as key_write:
            key_write.write(pem)
        with open(f"{key_out}.pubk", "wb") as key_write:
            key_write.write(pem2)

        return

    def read_pvk(self, pvk):
        with open(pvk, "rb") as key_file:
            self.private_key = serialization.load_pem_private_key(
                    key_file.read(),
                    password=None,
                    backend=default_backend()
                    )
        return

    def read_pubk(self, pubk):
        with open(pubk, "rb") as key_file:
            self.public_key = serialization.load_pem_public_key(
                    key_file.read(),
                    backend=default_backend()
                    )

    def make_fernet_key(self):
        self.encrypted_fernet = self.public_key.encrypt(
                Fernet.generate_key(),
                padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
                             algorithm=hashes.SHA256(),
                             label=None)
                )
        return

    def store_fernet_key(self, fernet_out):
        with open(fernet_out, "wb") as key_out:
            key_out.write(self.encrypted_fernet)
        return

    def read_fernet_key(self, fernet_file):
        with open(fernet_file, "rb") as key_file:
            self.encrypted_fernet = key_file.read()
        return

    def unlock_fernet(self):
        return self.private_key.decrypt(
                self.encrypted_fernet,
                padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
                             algorithm=hashes.SHA256(),
                             label=None)
                )

    def encrypt_with_fernet(self, unlocked):
        with open(unlocked, "rb") as infile:
            uncrypted = infile.read()
        encrypted = Fernet(self.unlock_fernet()).encrypt(uncrypted)
        with open(f"{unlocked}.locked", "wb") as outfile:
            outfile.write(encrypted)
        return

    def decrypt_with_fernet(self, locked):
        with open(locked, "rb") as infile:
            encrypted = infile.read()
        uncrypted = Fernet(self.unlock_fernet()).decrypt(encrypted)
        with open(f"{locked}.unlocked", "wb") as outfile:
            outfile.write(uncrypted)
        return


if __name__ == "__main__":
    import sys

    # ==Hashing==

    # Read pre-made hashtable.
    with open(sys.argv[2]) as hash_table:
        hash_values = hash_table.readlines()

    # Read corresponding varcon file.
    with open(sys.argv[1]) as varcon_in:
        varcons = varcon_in.readlines()

    # Parse text from each file.
    hash_values = [line.strip().split("\t") for line in hash_values if not line.startswith("#")]
    varcon_header = [line for line in varcons if line.startswith("#")]
    varcons = [line.split("\t") for line in varcons if not line.startswith("#")]

    outlines = ["\t".join([x[0], y[2], *x[2:]]) for x, y in zip(varcons, hash_values)]

    # Write hashed varcon file.
    with open(sys.argv[1] + ".hashed", "w") as hash_out:
        hash_out.writelines(varcon_header)
        hash_out.writelines(outlines)

    # ==Encryption==

    locker = Locker()

    locker.make_key_pair()
    locker.make_fernet_key()

    locker.store_key_pair("key_pair")
    locker.store_fernet_key("fernet.key")

    locker.encrypt_with_fernet(sys.argv[2])

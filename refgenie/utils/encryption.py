import json
import os
from base64 import b64decode, b64encode

from cryptography.fernet import Fernet

from refgenie.const import DEFAULT_ENCRYPTION_KEY
from refgenie.logger import logger


def get_encryption_key() -> bytes:
    """
    Get the encryption key from environment variable or use the default key.
    For production use, it's recommended to set REFGENIE_ENCRYPTION_KEY.

    Returns:
        bytes: The encryption key
    """
    if not (key := os.environ.get("REFGENIE_ENCRYPTION_KEY")):
        print(
            "WARNING: Using default encryption key for credential storage. "
            "For production use, set REFGENIE_ENCRYPTION_KEY. "
            'Generate one with: python -c "from cryptography.fernet import Fernet; '
            'print(Fernet.generate_key().decode())"'
        )
        key = DEFAULT_ENCRYPTION_KEY
    return b64decode(key)


def encrypt_dict(credentials: dict[str, str] | None) -> str | None:
    """
    Encrypt credentials dictionary using Fernet (symmetric encryption).

    Args:
        credentials: Dictionary of credentials to encrypt

    Returns:
        str | None: Encrypted credentials as a base64 string, or None if input is None
    """
    if credentials is None:
        return None

    if not credentials:
        return ""

    try:
        f = Fernet(get_encryption_key())
        json_data = json.dumps(credentials)
        encrypted_data = f.encrypt(json_data.encode())
        return b64encode(encrypted_data).decode("utf-8")
    except Exception as e:
        logger.error(f"Failed to encrypt credentials: {str(e)}")
        raise


def decrypt_credentials(encrypted_data: str | None) -> dict[str, str] | None:
    """
    Decrypt credentials string back to dictionary.

    Args:
        encrypted_data: Encrypted credentials string

    Returns:
        dict[str, str] | None: Decrypted credentials dictionary, or None if input is None
    """
    if encrypted_data is None:
        return None

    if not encrypted_data:
        return {}

    try:
        f = Fernet(get_encryption_key())
        decrypted_data = f.decrypt(b64decode(encrypted_data))
        return json.loads(decrypted_data)
    except Exception as e:
        logger.error(f"Failed to decrypt credentials: {str(e)}")
        raise

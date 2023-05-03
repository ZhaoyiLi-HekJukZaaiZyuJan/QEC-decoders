from itertools import chain
from typing import List

def hamming_encode(bits: List[int], r: int) -> List[int]:
    """Encode a list of bits using a Hamming code with r parity bits."""
    n = len(bits)
    k = n - r
    # Allocate space for encoded bits
    encoded_bits = [0] * (n + r)
    # Copy data bits into encoded bits, leaving space for parity bits
    j = 0
    for i in range(n + r):
        if i + 1 not in (1, 2, 4, 8, 16, 32, 64):  # Skip parity bit positions
            encoded_bits[i] = bits[j]
            j += 1
    # Compute parity bits
    for i in range(r):
        # Compute parity bit for each bit position that is a power of 2
        parity_bit = 0
        for j in range(1, n + r + 1):
            if j & (2**i):
                parity_bit ^= encoded_bits[-j]
        # Set parity bit in encoded bits
        encoded_bits[-(2**i)] = parity_bit
    return encoded_bits

def hamming_decode(bits: List[int], r: int) -> List[int]:
    """Decode a list of bits encoded using a Hamming code with r parity bits."""
    n = len(bits) - r
    k = n - r
    # Allocate space for decoded bits
    decoded_bits = [0] * n
    # Compute syndrome
    syndrome = 0
    for i in range(r):
        # Compute parity bit for each bit position that is a power of 2
        parity_bit = 0
        for j in range(1, n + r + 1):
            if j & (2**i):
                parity_bit ^= bits[-j]
        # Add parity bit to syndrome
        syndrome |= (parity_bit << i)
    # If syndrome is non-zero, correct error
    if syndrome:
        # Flip bit at position corresponding to syndrome
        bits[-(syndrome)] ^= 1
    # Copy data bits into decoded bits
    j = 0
    for i in range(n + r):
        if i + 1 not in (1, 2, 4, 8, 16, 32, 64):  # Skip parity bit positions
            decoded_bits[j] = bits[i]
            j += 1
    return decoded_bits

def hamming_concatenated_encode(bits: List[int], r: int, n_blocks: int) -> List[int]:
    """Encode a list of bits using a concatenated Hamming code with n_blocks blocks of r parity bits."""
    block_size = len(bits) // n_blocks
    block_encoded_bits = []
    for i in range(n_blocks):
        block_bits = bits[i*block_size : (i+1)*block_size]
        block_encoded_bits += hamming_encode(block_bits, r)
    return block_encoded_bits

def hamming_concatenated_decode(bits: List[int], r: int, n_blocks: int) -> List[int]:
    """Decode a list of bits encoded using a concatenated Hamming code with n_blocks blocks of r parity bits."""
    block_size = len(bits) // n_blocks
    block_decoded_bits = []
    for i in range(n_blocks):
        block_bits = bits[i*block_size : (i+1)*block_size]
        block_decoded_bits += hamming
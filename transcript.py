import hashlib
import struct

class Transcript:
    """
    Handles the protocol transcript and deterministic challenge generation for Fiat-Shamir transform.
    
    This class maintains the state of all messages exchanged during the protocol
    and generates deterministic challenges that depend on the transcript history.
    """
    
    def __init__(self, label, F):
        """
        Initialize a new transcript with a label.
        
        Args:
            label: A string label for this transcript instance
            F: The finite field used in the protocol
        """
        self.label = label
        self.F = F
        self.state = hashlib.sha256(label.encode()).digest()
    
    def append_message(self, message_label, message_data):
        """
        Append a message to the transcript.
        
        Args:
            message_label: A string label for this message
            message_data: The message data to append (will be serialized)
        """
        
        # Update the state
        self._update_state(message_label, self._serialize(message_data))
    
    def get_challenge(self, label):
        """
        Generate a deterministic challenge field element based on the current transcript.
        
        Args:
            label: A label for this challenge
            
        Returns:
            A field element derived from the current transcript state
        """
        # Update state with the challenge label
        challenge_state = hashlib.sha256(self.state + label.encode()).digest()
        
        # Convert the first 32 bytes of hash to an integer mod field order
        challenge_int = int.from_bytes(challenge_state, byteorder='big')
        challenge = self.F(challenge_int)
        
        # Update the transcript state with this challenge
        self._update_state(label, challenge_state)
        
        return challenge
    
    def _serialize(self, data):
        """
        Serialize various data types for inclusion in the transcript.
        
        Args:
            data: The data to serialize
            
        Returns:
            Bytes representation of the data
        """
        if isinstance(data, str):
            return data.encode()
        elif isinstance(data, int):
            return struct.pack(">q", data)  # 8-byte big-endian integer
        elif isinstance(data, bytes):
            return data
        elif isinstance(data, list):
            # For lists, serialize each element and concatenate
            result = b''
            for item in data:
                result += self._serialize(item)
            return result
        elif hasattr(data, '_repr_'):
            # For SageMath objects with string representation
            return str(data).encode()
        else:
            # Default: convert to string and encode
            return str(data).encode()
    
    def _update_state(self, label, data):
        """
        Update the internal state with a new message.
        
        Args:
            label: Label for the message
            data: Serialized message data
        """
        # Update state: H(state || label || data)
        hasher = hashlib.sha256()
        hasher.update(self.state)
        hasher.update(label.encode())
        hasher.update(data)
        self.state = hasher.digest()
package blazmass.dbindex;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

/**
 * Dynamic byte buffer utility
 *
 * @author Adam
 */
public class DynByteBuffer {

    private static int DEFAULT_CAPACITY = 4 * 4; //1 sequence = 4 * 4 bytes
    private int curSize = 0;
    private int curCapacity = 16;
    private byte[] data;

    DynByteBuffer() {
        data = new byte[curCapacity];
    }

    void add(byte[] toAdd) {
        final int toAddLen = toAdd.length;
        if (curCapacity < curSize + toAddLen) {
            int newCapacity = (curCapacity + toAddLen) * 3;
            byte[] newData = new byte[newCapacity];
            System.arraycopy(data, 0, newData, 0, curSize);
            curCapacity = newCapacity;
            data = newData;
        }

        //insert new data
        System.arraycopy(toAdd, 0, data, curSize, toAddLen);
        curSize += toAddLen;
    }

    void clear() {
        curSize = 0;
        curCapacity = DEFAULT_CAPACITY * 4;
        data = new byte[curCapacity];
        
    }
    int getSize() {
        return curSize;
    }

    byte[] getData() {
        return Arrays.copyOfRange(data, 0, curSize);
    }

    public static byte[] toByteArray(int myInteger) {
        return ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN).putInt(myInteger).array();
    }

    public static int toInt(byte[] byteBarray) {
        return ByteBuffer.wrap(byteBarray).order(ByteOrder.LITTLE_ENDIAN).getInt();
    }
    
    public static byte[] toByteArray(float theFloat) {
        return ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN).putFloat(theFloat).array();
    }

    public static float toFloat(byte[] byteBarray) {
        return ByteBuffer.wrap(byteBarray).order(ByteOrder.LITTLE_ENDIAN).getFloat();
    }
}

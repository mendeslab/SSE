package SSE;

/** A numerical kernel failure, not a zero-probability model state. No retry is performed. */
public class QuaSSEKernelException extends ArithmeticException {
    public QuaSSEKernelException(String message) {
        super(message);
    }

    public QuaSSEKernelException(String message, QuaSSEKernelException cause) {
        super(message);
        initCause(cause);
    }
}

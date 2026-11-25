package SSE;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.ArrayDeque;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.function.BiFunction;

/** Schedules ordinary branches bottom-up. The caller finishes the root after run returns. */
final class QuaSSEBranchTraversal implements AutoCloseable {
    private final ThreadPoolExecutor executor;
    private final int workers;

    // Idle threads may expire if a host application discards a likelihood without closing it.
    QuaSSEBranchTraversal(int workers) {
        this.workers = workers;
        executor = new ThreadPoolExecutor(workers, workers, 30, TimeUnit.SECONDS, new LinkedBlockingQueue<>(), r -> {
            Thread t = new Thread(r, "QuaSSE branch worker");
            t.setDaemon(true);
            return t;
        });
        executor.allowCoreThreadTimeOut(true);
    }

    // The coordinator owns the FIFO, counts completed children and leases a workspace slot to
    // each active branch. Future.get publishes both child rows before their parent is prepared.
    void run(Tree tree, BiFunction<Node, Integer, Runnable> prepare) {
        var ready = new ArrayDeque<>(tree.getExternalNodes());
        var free = new ArrayDeque<Integer>(workers);
        for (int i = 0; i < workers; i++) free.add(i);
        var completed = new ExecutorCompletionService<Integer>(executor);
        Node[] active = new Node[workers];
        int[] childrenDone = new int[tree.getNodeCount()];
        int pending = 0;
        boolean interrupted = false;
        try {
            while (!ready.isEmpty() || pending > 0) {
                while (!ready.isEmpty() && !free.isEmpty()) {
                    int slot = free.remove();
                    Node node = ready.remove();
                    Runnable integrate = prepare.apply(node, slot);
                    active[slot] = node;
                    completed.submit(integrate, slot);
                    pending++;
                }
                var result = completed.take();
                pending--;
                int slot = result.get();
                Node parent = active[slot].getParent();
                if (!parent.isRoot() && ++childrenDone[parent.getNr()] == 2) ready.add(parent);
                free.add(slot);
            }
        } catch (InterruptedException e) {
            interrupted = true;
            throw new IllegalStateException("Interrupted during QuaSSE branch integration", e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof Error error) throw error;
            if (cause instanceof RuntimeException error) throw error;
            throw new IllegalStateException(cause);
        } finally {
            // Even an interrupted native call can still write into the proposal's rows. Wait for
            // every dispatched task to finish before the caller invalidates keys or restores rows.
            while (pending > 0) {
                try {
                    completed.take();
                    pending--;
                } catch (InterruptedException e) {
                    interrupted = true;
                }
            }
            if (interrupted) Thread.currentThread().interrupt();
        }
    }

    @Override
    public void close() {
        executor.shutdown();
    }
}

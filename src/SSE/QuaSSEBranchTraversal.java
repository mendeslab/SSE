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
import java.util.function.IntConsumer;

/** Schedules ordinary branches bottom-up. The caller finishes the root after run returns. */
final class QuaSSEBranchTraversal implements AutoCloseable {
    private final ThreadPoolExecutor executor;
    private final int workers;

    // Idle threads may expire if a host application discards a likelihood without closing it.
    QuaSSEBranchTraversal(int workers) {
        this.workers = workers;
        if (workers == 1) {
            executor = null;
            return;
        }
        executor = new ThreadPoolExecutor(workers, workers, 30, TimeUnit.SECONDS, new LinkedBlockingQueue<>(), r -> {
            Thread t = new Thread(r, "QuaSSE branch worker");
            t.setDaemon(true);
            return t;
        });
        executor.allowCoreThreadTimeOut(true);
    }

    // Visit children in completion order, preserving tree order for the actual merge. A hit
    // (prepare returns null) counts as a finished child. A recomputed child also marks the parent
    // for integration. The caller owns the queue and leases each worker its slot until get().
    void run(Tree tree, BiFunction<Node, Boolean, IntConsumer> prepare) {
        var ready = new ArrayDeque<>(tree.getExternalNodes());
        var free = new ArrayDeque<Integer>(workers);
        for (int i = 0; i < workers; i++) free.add(i);
        var completed = executor == null ? null : new ExecutorCompletionService<Integer>(executor);
        Node[] active = new Node[workers];
        int[] childrenDone = new int[tree.getNodeCount()];
        boolean[] childChanged = new boolean[tree.getNodeCount()];
        int pending = 0;
        boolean interrupted = false;
        try {
            while (!ready.isEmpty() || pending > 0) {
                Node node;
                boolean changed;
                if (!ready.isEmpty() && !free.isEmpty()) {
                    node = ready.remove();
                    IntConsumer integrate = prepare.apply(node, childChanged[node.getNr()]);
                    changed = integrate != null;
                    if (changed && completed != null) {
                        int slot = free.remove();
                        active[slot] = node;
                        completed.submit(() -> integrate.accept(slot), slot);
                        pending++;
                        continue;
                    }
                    if (changed) integrate.accept(0);
                } else {
                    // The ready queue is empty or all slots are occupied. In either case
                    // at least one task is outstanding, since the loop ends when both are empty.
                    var result = completed.take();
                    pending--;
                    int slot = result.get();
                    node = active[slot];
                    changed = true;
                    free.add(slot);
                }
                Node parent = node.getParent();
                if (!parent.isRoot()) {
                    int i = parent.getNr();
                    childChanged[i] |= changed;
                    if (++childrenDone[i] == 2) ready.add(parent);
                }
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
        if (executor != null) executor.shutdown();
    }
}

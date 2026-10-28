"""Minimal rooted tree with last-common-ancestor queries.

Deliberately independent of ete3: the taxonomy backend is one optional loader,
so the core algorithms can be tested without a 2.9M-node NCBI dump. A tree is
just a parent map plus a leaf set.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence


class Tree:
    """A rooted tree over hashable node ids."""

    def __init__(self, parent: Mapping[str, str | None]):
        self.parent: dict[str, str | None] = dict(parent)
        roots = [n for n, p in self.parent.items() if p is None]
        if len(roots) != 1:
            raise ValueError(f"expected exactly one root, found {len(roots)}: {roots[:5]}")
        self.root = roots[0]

        self.children: dict[str, list[str]] = {n: [] for n in self.parent}
        for node, par in self.parent.items():
            if par is not None:
                if par not in self.children:
                    raise ValueError(f"parent {par!r} of {node!r} is not a node")
                self.children[par].append(node)

        self.leaves: list[str] = sorted(n for n, kids in self.children.items() if not kids)
        self._depth = self._compute_depths()
        self._descendant_leaves = self._compute_descendant_leaves()

    def _compute_depths(self) -> dict[str, int]:
        depth: dict[str, int] = {}
        for node in self.parent:
            chain = []
            cur: str | None = node
            while cur is not None and cur not in depth:
                chain.append(cur)
                cur = self.parent[cur]
            base = 0 if cur is None else depth[cur] + 1
            for offset, name in enumerate(reversed(chain)):
                depth[name] = base + offset
        return depth

    def _compute_descendant_leaves(self) -> dict[str, frozenset[str]]:
        out: dict[str, frozenset[str]] = {}

        def visit(node: str) -> frozenset[str]:
            if node in out:
                return out[node]
            kids = self.children[node]
            acc = frozenset([node]) if not kids else frozenset().union(*(visit(k) for k in kids))
            out[node] = acc
            return acc

        # iterative post-order to stay safe on deep taxonomies
        stack = [(self.root, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                kids = self.children[node]
                out[node] = frozenset([node]) if not kids else frozenset().union(*(out[k] for k in kids))
            else:
                stack.append((node, True))
                stack.extend((k, False) for k in self.children[node])
        return out

    def depth(self, node: str) -> int:
        return self._depth[node]

    def descendant_leaves(self, node: str) -> frozenset[str]:
        return self._descendant_leaves[node]

    def ancestors(self, node: str) -> list[str]:
        """`node` first, then each ancestor up to and including the root."""
        chain, cur = [], node
        while cur is not None:
            chain.append(cur)
            cur = self.parent[cur]
        return chain

    def lca(self, nodes: Iterable[str]) -> str:
        nodes = list(nodes)
        if not nodes:
            raise ValueError("lca() of an empty set is undefined")
        missing = [n for n in nodes if n not in self.parent]
        if missing:
            raise KeyError(f"unknown node(s): {missing[:5]}")
        current = nodes[0]
        for other in nodes[1:]:
            current = self._lca_pair(current, other)
        return current

    def _lca_pair(self, a: str, b: str) -> str:
        da, db = self._depth[a], self._depth[b]
        while da > db:
            a = self.parent[a]  # type: ignore[assignment]
            da -= 1
        while db > da:
            b = self.parent[b]  # type: ignore[assignment]
            db -= 1
        while a != b:
            a = self.parent[a]  # type: ignore[assignment]
            b = self.parent[b]  # type: ignore[assignment]
        return a

    def nodes(self) -> list[str]:
        return list(self.parent)

    @classmethod
    def from_edges(cls, edges: Sequence[tuple[str, str]], root: str) -> "Tree":
        """Build from (parent, child) pairs."""
        parent: dict[str, str | None] = {root: None}
        for par, child in edges:
            parent[child] = par
            parent.setdefault(par, None) if par == root else parent.setdefault(par, parent.get(par))
        for par, _child in edges:
            if par not in parent:
                raise ValueError(f"edge parent {par!r} never appears as a child and is not the root")
        return cls(parent)


def flat_two_level_tree(groups: Mapping[str, Sequence[str]], root: str = "root") -> Tree:
    """root -> group -> leaves. Handy for tests and simulation."""
    parent: dict[str, str | None] = {root: None}
    for group, members in groups.items():
        parent[group] = root
        for leaf in members:
            parent[leaf] = group
    return Tree(parent)

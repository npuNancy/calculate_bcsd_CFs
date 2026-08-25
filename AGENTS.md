# Codex Instructions

## Communication

- 默认用英文思考，用中文回答，除非我明确要求用英文回答。
- 回答要直接、具体，避免空泛解释。
- 修改代码前先理解现有结构和约定。
- 每次回复我时，都称呼我为 `小凯`，并且在回复的最后加上 `希望对你有帮助，小凯！`。

## Engineering Defaults

- 优先遵循项目已有风格，不引入不必要抽象。
- 修改后尽量运行相关测试、lint 或类型检查。
- 不要重构无关代码。
- 遇到脏工作区时，不要回滚我已有的改动。
- 当前工作区的 python 使用 uv 环境: `.venv`
- `RQ3/` 的主脚本和 `RQ3/tmp/` 共享数据入口当前包含临时代码：因上游风电数据被错误放大 10 倍，脚本会将相关风电发电量与损失能量字段统一乘以 0.1；上游数据修复后需移除并重跑结果。

### Scope Discipline

Implement the requested change, not the story behind the change.

Do:
- Make the smallest complete change.
- Keep names, APIs, and docs focused on the final desired state.

Don't:
- Add features, rules, comments, or documentation explaining rejected ideas.
- Preserve removed concepts in names (e.g., "no-X", "without-X").
- Convert a correction into a general design principle.

A fix is a fix, not a new product feature.

## Safety

- 不要执行破坏性 git 命令，除非我明确要求。
- 新增依赖前说明原因。
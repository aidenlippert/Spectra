# Implementation

- Build the smallest complete behavior that satisfies the request, including its required technical capability and refusal paths.
- Reuse or remove existing code before adding another implementation. Prefer the standard library, existing dependencies, and the nearest suitable module over new frameworks, services, registries, wrappers, or configuration.
- Add an abstraction for a current external boundary, a replaceable implementation, or multiple real callers. Keep interfaces narrow, modules deep, data flow direct, and related behavior local.
- Put working behavior on the production path with a named caller. Schemas, mocks, placeholders, and speculative infrastructure do not satisfy a request for working product behavior.
